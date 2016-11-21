
#include <mpi.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <cstring>

#include "Reaction.hh"
#include "VectorDouble32.hh"
#include "reactionFactory.hh"
#include "Anatomy.hh"
#include "ThreadServer.hh"
#include "object.h"
#include "object_cc.hh"

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


int main(int argc, char* argv[])
{

   int npes, mype;
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &npes);
   MPI_Comm_rank(MPI_COMM_WORLD, &mype);

   struct gengetopt_args_info params;
   cmdline_parser(argc, argv, &params);

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

   //read in the object file
   for (int iobject=0; iobject<params.object_given; iobject++) 
   {
      object_compilefile(params.object_arg[iobject]);
   }
   
   //figure out the name of the object
   string objectName;
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

   //create the dependent variables.
   Anatomy anatomy;
   {
      anatomy.setGridSize(nCells,1,1);
      vector<AnatomyCell>& cells = anatomy.cellArray();
      cells.reserve(nCells);
      for (int ii=0; ii<nCells; ii++)
      {
         AnatomyCell tmp;
         tmp.gid_ = ii;
         tmp.cellType_ = params.anatomy_tag_arg;
         cells.push_back(tmp);
      }
   }
   ThreadServer& threadServer = ThreadServer::getInstance();
   ThreadTeam threads = threadServer.getThreadTeam(vector<unsigned>());
   Reaction* reaction = reactionFactory(objectName,dt,
                                        anatomy, threads,
                                        vector<string>());

   //initialize the ionic model
   VectorDouble32 Vm(nCells);
   vector<double> iStim(nCells);
   VectorDouble32 dVm(nCells);
   initializeMembraneState(reaction, objectName, Vm);

   //prepare for checkpoint output and extra columns
   vector<string> fieldNames;
   vector<string> fieldUnits;
   reaction->getCheckpointInfo(fieldNames,fieldUnits);
   vector<int> fieldHandles(fieldNames.size());
   for (int istate=0; istate<fieldNames.size(); ++istate)
   {
      fieldHandles[istate] = reaction->getVarHandle(fieldNames[istate]);
   }

   vector<int> extraColumns(params.add_column_given);
   for (int iextra=0; iextra<extraColumns.size(); ++iextra) 
   {
      //find the state name
      int foundState = -1;
      for (int istate=0; istate<fieldNames.size(); ++istate) 
      {
         if (fieldNames[istate] == params.add_column_arg[iextra]) 
         {
            foundState = istate;
            break;
         }
      }
      if (foundState == -1)
      {
         fprintf(stderr, "Couldn't find state '%s' in the model\n", params.add_column_arg[iextra]);
         return __LINE__;
      }
      extraColumns[iextra] = fieldHandles[foundState];
   }
   
   //using a forever loop here so we can repeat outputs on the last iteration.
   //otherwise, we'd use a for loop here.
   int itime=0;
   while(1)
   {
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
            fprintf(file, "%s REACTION { initialState = %s; }\n",
                    objectName.c_str(), objectName.c_str());
            fprintf(file, "%s SINGLECELL {\n", objectName.c_str());
            fprintf(file, "  method = %s;\n", reaction->methodName().c_str());
            fprintf(file, "  Vm = %.17g mV;\n", Vm[0]);
            for (int istate=0; istate<fieldNames.size(); ++istate)
            {
               fprintf(file, "  %s = %.17g %s;\n",
                       fieldNames[istate].c_str(),
                       reaction->getValue(0,fieldHandles[istate]),
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
         printf("%21.17g %21.17g", timeline.realTimeFromTimestep(itime), Vm[0]);
         for (int iextra=0; iextra<extraColumns.size(); ++iextra) 
         {
            printf(" %21.17g", reaction->getValue(0,extraColumns[iextra]));
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
      for (int ii=0; ii<nCells; ii++)
      {
         iStim[ii]=0;
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
            iStim[ii] = -params.stim_strength_arg;
            timestepsLeftInStimulus--;
         }
      }

      //calc() dVm and update the internal states
      reaction->calc(dt, Vm, iStim, dVm);

      //update Vm and apply stimulus
      for (int ii=0; ii<nCells; ii++)
      {
         //use a negative sign here to undo the negative we had above.
         Vm[ii] += (dVm[ii] - iStim[ii]) * dt;
      }

      itime++;
   }


   MPI_Finalize();
   return 0;
}
