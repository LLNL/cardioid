
#include <mpi.h>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>

#include "Reaction.hh"
#include "VectorDouble32.hh"
#include "reactionFactory.hh"
#include "Anatomy.hh"
#include "ThreadServer.hh"
#include "object.h"
#include "object_cc.hh"
#include "DeviceFor.hh"

#include "singleCellOptions.h"

using namespace std;

MPI_Comm COMM_LOCAL = MPI_COMM_WORLD;

/*!
  Keeps track of real time versus simulation ticks.
 */
class Timeline
{
 public:
   Timeline(double dt, double duration)
   {
      maxTimesteps_ = round(duration/dt);
      dt_ = duration/maxTimesteps_;
   }
   int maxTimesteps() const { return maxTimesteps_; };
   double dt() const { return dt_; }
   double maxTime() const { return dt_*maxTimesteps_; }
   double realTimeFromTimestep(int timestep) const
   {
      return timestep*dt_;
   }
   int timestepFromRealTime(double realTime) const
   {
      return round(realTime/dt_);
   }

 private:
   double dt_;
   int maxTimesteps_;
};

class FunctionFromFile
{
 public:
   double operator()(const double x) const
   {
      if (valuesFromTime_.empty()) { assert(0 && "Need to populate the map!" ); }
      std::map<double,double>::const_iterator upper=valuesFromTime_.upper_bound(x);
      if (upper == valuesFromTime_.begin()) { return upper->second; }
      std::map<double,double>::const_iterator lower=upper;
      lower--;
      if (upper == valuesFromTime_.end()) { return lower->second; }
      if (lower->first == upper->first) { return (lower->second+upper->second)/2; }
      return ((x-lower->first)*upper->second + (upper->first-x)*lower->second)/(upper->first-lower->first);
   }
   void addData(const double x, const double y)
   {
      valuesFromTime_.insert(make_pair(x,y));
   }
   
 private:
   std::multimap<double,double> valuesFromTime_;
};


int main(int argc, char* argv[])
{

   int npes, mype;
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &npes);
   MPI_Comm_rank(MPI_COMM_WORLD, &mype);

#ifdef USE_CUDA
   // A ugly way to trigger the default CUDA context
   int *d_i;
   cudaMalloc(&d_i, sizeof(int));
   cudaFree(d_i);
#endif

   struct gengetopt_args_info params;
   cmdline_parser(argc, argv, &params);
   if (!params.object_given && !params.method_given)
   {
      cmdline_parser_print_help();
      cerr << "ERROR: Must specify either at least one object or a method name.\n";
      return -1;
   }
   else if (params.object_given && params.method_given)
   {
      cmdline_parser_print_help();
      cerr << "ERROR: Cannot specify both a method name and an object file-- pick one.\n";
      return -1;
   }
   
   double dt = params.dt_arg;

   double duration;
   if (params.duration_given)
   {
      duration = params.duration_arg;
   }
   else
   {
      duration = params.s1_offset_arg + params.s1_count_arg*params.s1_bcl_arg;
   }

   Timeline timeline(dt, duration);


   //Stimulus setup
   vector<int> stimTimesteps;
   if (params.stim_at_given)
   {
      for(int istim=0; istim<params.stim_at_given; istim++)
      {
         stimTimesteps.push_back(timeline.timestepFromRealTime(params.stim_at_arg[istim]));
      }
   }
   bool s1_set_on_command_line = params.s1_count_given || params.s1_bcl_given || params.s1_offset_given;
   if (s1_set_on_command_line || !params.stim_at_given)
   {
      for (int ibeat=0; ibeat<params.s1_count_arg; ibeat++)
      {
         stimTimesteps.push_back(timeline.timestepFromRealTime
                                 (params.s1_offset_arg + ibeat*params.s1_bcl_arg));
      }
   }
   sort(stimTimesteps.begin(), stimTimesteps.end());

   int timestepsInEachStimulus=timeline.timestepFromRealTime(params.stim_duration_arg);
   double stimulusStrength = params.stim_strength_arg;
   int timestepsLeftInStimulus=0;
   vector<int>::iterator nextStimulus=stimTimesteps.begin();


   //Output setup
   int outputTimestepInterval = timeline.timestepFromRealTime(params.output_dt_arg);

   bool shouldWriteState=false;
   int writeStateTimestep=timeline.maxTimesteps();
   if (params.save_state_time_given || params.save_state_file_given)
   {
      shouldWriteState=true;
      if (params.save_state_time_given)
      {
         writeStateTimestep=timeline.timestepFromRealTime(params.save_state_time_arg);
      }
   }
   
   //create the ionic model
   const int nCells = params.num_points_arg;

   string objectName;
   if (params.method_given)
   {
      string objectFile = string("default REACTION { method=")+params.method_arg+"; }";
      object_compilestring(const_cast<char*>(objectFile.c_str()));
      objectName = "default";
   }
   else
   {
      //read in the object file
      for (int iobject=0; iobject<params.object_given; iobject++) 
      {
         object_compilefile(params.object_arg[iobject]);
      }
   
      if (params.reaction_name_given) 
      {
         objectName = params.reaction_name_arg;
         objectFind(objectName, "REACTION");
      }
      else
      {
         /* 
            Unlike most uses of the object database, we'd like users to be
            able to name the top level object here whatever they'd like, so
            they can use the same reaction object in actual simulation runs
         */
         string objectClass;
         OBJECTFILE file;
         file = object_fopen(params.object_arg[0], "r");
         if (file.file == NULL)
         {
            fprintf(stderr, "Object file could not be opened for reading.\n");
            return __LINE__;
         }
         char* g_line = object_read(file);
         if (g_line == NULL)
         {
            fprintf(stderr, "Object file is empty.\n");
            return __LINE__;
         }
         char* line = strdup(g_line);
         OBJECT obj;
         int rc = object_lineparse(line, &obj);
         if (rc >= 2)
         {
            fprintf(stderr, "Other object parsing error?\n");
            return __LINE__;
         }
         objectName = obj.name;
         objectClass = obj.objclass;

         free(obj.name);
         free(obj.objclass);
         free(obj.value);
         free(line);
         object_fclose(file);
         if (objectClass != "REACTION")
         {
            fprintf(stderr, "First object in file must be a REACTION!\n");
            return __LINE__;
         }
      }
   }

   //create the dependent variables.
   ThreadServer& threadServer = ThreadServer::getInstance();
   ThreadTeam threads = threadServer.getThreadTeam(vector<unsigned>());
   Reaction* reaction = reactionFactory(objectName,dt,
                                        nCells, threads);

   //initialize the ionic model
   lazy_array<double> VmTransport;
   VmTransport.resize(nCells);
   lazy_array<double> iStimTransport;
   iStimTransport.resize(nCells);
   lazy_array<double> dVmTransport;
   dVmTransport.resize(nCells);
   lazy_array<int> indexTransport;
   indexTransport.resize(nCells);
   {
      wo_array_ptr<int> indexArray = indexTransport.useOn(CPU);
      for (int ii=0; ii<nCells; ii++)
      {
         indexArray[ii] = ii;
      }
   }
   
   initializeMembraneState(reaction, objectName, indexTransport, VmTransport);

   //prepare for checkpoint output and extra columns
   vector<string> fieldNames;
   vector<string> fieldUnits;
   reaction->getCheckpointInfo(fieldNames,fieldUnits);
   vector<int> fieldHandles(fieldNames.size());
   for (int istate=0; istate<fieldNames.size(); ++istate)
   {
      fieldHandles[istate] = reaction->getVarHandle(fieldNames[istate]);
   }

   vector<int> extraColumns;
   vector<string> extraColumnHeaders;
   if (params.add_all_state_given)
   {
      for (int istate=0; istate<fieldNames.size(); ++istate) 
      {
         extraColumns.push_back(fieldHandles[istate]);
         extraColumnHeaders.push_back(fieldNames[istate]);
      }
   }
   bool allColumnsAdded = true;
   for (int iextra=0; iextra<params.add_column_given; ++iextra)
   {
      int handle = reaction->getVarHandle(params.add_column_arg[iextra]);
      if (handle != -1)
      {
         extraColumns.push_back(handle);
         extraColumnHeaders.push_back(params.add_column_arg[iextra]);
      }
      else
      {
         fprintf(stderr, "Can't add column '%s': can't find the handle.\n",params.add_column_arg[iextra]);
         allColumnsAdded = false;
      }
   }
   if (!allColumnsAdded)
   {
      exit(1);
   }

   if (params.add_header_given)
   {
      printf("                     time                        Vm");
      for (int ii=0; ii<extraColumnHeaders.size(); ++ii)
      {
         printf(" %25s",extraColumnHeaders[ii].c_str());
      }
      printf("\n");
   }

   const int VmKey = -99;
   vector<int> clampHandles;
   vector<FunctionFromFile> clampFunctions;
   if (params.clamp_file_given)
   {
      ifstream clampFile(params.clamp_file_arg);
      string temp;
      getline(clampFile,temp);
      stringstream headerLine(temp);
      string timeName;
      headerLine >> timeName;
      assert(timeName == "time" || timeName == "t");
      while (1)
      {
         string columnName;
         headerLine >> columnName;
         if (!headerLine) { break; }
         int handle;
         if (columnName == "Vm")
         {
            handle = VmKey;
         }
         else
         {
            handle = reaction->getVarHandle(columnName);
            assert(handle >= 0 && "Invalid column name!");
         }
         clampHandles.push_back(handle);
      }
      vector<double> timestepValues(clampHandles.size());
      clampFunctions.resize(clampHandles.size());
      while (1)
      {
         string temp;
         getline(clampFile,temp);
         stringstream thisLine(temp);
         double time;
         thisLine >> time;
         for (int iclamp=0; iclamp<clampHandles.size(); iclamp++)
         {
            thisLine >> timestepValues[iclamp];
         }
         if (!thisLine) { break; }
         for (int iclamp=0; iclamp<clampHandles.size(); iclamp++)
         {
            clampFunctions[iclamp].addData(time,timestepValues[iclamp]);
         }
      }
   }


   //using a forever loop here so we can repeat outputs on the last iteration.
   //otherwise, we'd use a for loop here.
   int itime=0;
   while(1)
   {
      //if we should clamp, do so.
      for (int iclamp=0; iclamp<clampHandles.size(); iclamp++)
      {
         int handle = clampHandles[iclamp];
         double value = clampFunctions[iclamp](timeline.realTimeFromTimestep(itime));
         for (int ii=0; ii<nCells; ii++)
         {
            if (handle == VmKey)
            {
               wo_array_ptr<double> Vm = VmTransport.useOn(CPU);
               Vm[ii] = value;
            }
            else
            {
               reaction->setValue(ii,handle,value);
            }
         }
      }
      
      //if we should checkpoint, do so.
      if (shouldWriteState && itime == writeStateTimestep)
      {
         FILE* file = fopen(params.save_state_file_arg, "w");
         if (file == NULL)
         {
            perror(("Can't open "+string(params.save_state_file_arg)+" for writing: ").c_str());
         }
         else
         {
            ro_array_ptr<double> Vm = VmTransport.useOn(CPU);
            fprintf(file, "%s REACTION { initialState = %s; }\n",
                    objectName.c_str(), objectName.c_str());
            fprintf(file, "%s SINGLECELL {\n", objectName.c_str());
            fprintf(file, "  method = %s;\n", reaction->methodName().c_str());
            fprintf(file, "  Vm = %.17g mV;\n", Vm[0]);
            for (int istate=0; istate<fieldNames.size(); ++istate)
            {
               fprintf(file, "  %s = %.17g %s;\n",
                       fieldNames[istate].c_str(),
                       reaction->getValue(0,fieldHandles[istate], Vm[0]),
                       fieldUnits[istate].c_str()
                      );
            }
            fprintf(file, "}\n");
            fclose(file);
         }
      }

      //if we should do output, do so.
      if ((itime % outputTimestepInterval) == 0)
      {
         //doIO();
         ro_array_ptr<double> Vm = VmTransport.useOn(CPU);
         printf("%21.17g %21.17g", timeline.realTimeFromTimestep(itime), Vm[0]);
         for (int iextra=0; iextra<extraColumns.size(); ++iextra) 
         {
            printf(" %21.17g", reaction->getValue(0,extraColumns[iextra],Vm[0]));
         }
         printf("\n");
      }

      //if we're done with the loop, exit.
      if (itime == timeline.maxTimesteps()) { break; }

      //any stimulii left?
      if (nextStimulus != stimTimesteps.end())
      {
         //if so, wait for it.
         if (itime == *nextStimulus)
         {
            ++nextStimulus;
            timestepsLeftInStimulus=timestepsInEachStimulus;
         }
      }
      {
         double stimAmount = 0;
         if (timestepsLeftInStimulus > 0)
         {
            /* Check the negative sign here.  Look at:
               
               startTimer(stimulusTimer);
               // add stimulus to dVmDiffusion
               for (unsigned ii=0; ii<sim.stimulus_.size(); ++ii)
                  sim.stimulus_[ii]->stim(sim.time_, dVmDiffusion);
               for (unsigned ii=0; ii<nLocal; ++ii)
                  iStim[ii] = -(dVmDiffusion[ii]);
               stopTimer(stimulusTimer);

               Notice that the stimulii all add their current to
               stimulus. If you assume diffusion is zero, then istim ends
               up negative.

            */
            stimAmount = -params.stim_strength_arg;
            timestepsLeftInStimulus--;
         }
         wo_array_ptr<double> iStim = iStimTransport.writeonly(DEFAULT_COMPUTE_SPACE);
         DEVICE_PARALLEL_FORALL(iStim.size(), ii,
                                iStim[ii] = -stimAmount);

      }

      {         
         //calc() dVm and update the internal states
         if (params.alternate_update_flag)
         {
            reaction->updateNonGate(dt, indexTransport, VmTransport, dVmTransport);
            reaction->updateGate(dt, indexTransport, VmTransport);
         }
         else
         {
            reaction->calc(dt, indexTransport, VmTransport, iStimTransport, dVmTransport);
         }
      }
      {
         //update Vm and apply stimulus
         rw_array_ptr<double> Vm = VmTransport.readwrite(DEFAULT_COMPUTE_SPACE);
         ro_array_ptr<double> iStim = iStimTransport.readonly(DEFAULT_COMPUTE_SPACE);
         ro_array_ptr<double> dVm = dVmTransport.readonly(DEFAULT_COMPUTE_SPACE);
         DEVICE_PARALLEL_FORALL(VmTransport.size(), ii,
                                Vm[ii] += dt*(dVm[ii]+iStim[ii]));
      }
      itime++;
   }


   MPI_Finalize();
   return 0;
}
