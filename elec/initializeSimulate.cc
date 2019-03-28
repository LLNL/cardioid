#include "initializeSimulate.hh"

#include <iostream>
#include <sstream>
#include <map>

#include "object_cc.hh"
#include "Simulate.hh"

#include "initializeAnatomy.hh"
#include "assignCellsToTasks.hh"
#include "diffusionFactory.hh"
#include "ReactionManager.hh"
#include "stimulusFactory.hh"
#include "Stimulus.hh"
#include "sensorFactory.hh"
#include "getRemoteCells.hh"
#include "Anatomy.hh"
#include "mpiUtils.h"
#include "PerformanceTimers.hh"
#include "hardwareInfo.h"
#include "pio.h"
#include "heap.h"
#include "LoadLevel.hh"

using namespace std;

namespace
{
   void checkForObsoleteKeywords(OBJECT* obj);
   void buildCoreList(unsigned& nCores, vector<unsigned>& cores);
   Long64 findGlobalMinGid(const Anatomy& anatomy);
   void writeTorusMap(MPI_Comm comm, const string& filename);
}



/*!
   @page obj_SIMULATE SIMULATE object

   The SIMULATE object is the master object of the simulation.  It
   specifies simulation parameters such as the time step, maximum loop
   count, etc, as well as the names of the other objects that define the
   reaction and diffusion models, and the stimulus and sensor
   protocols.
   
   @beginkeywords
   @kw{anatomy, The name of the ANATOMY object for this simulation., anatomy}
   @kw{checkpointRate, The rate (in time steps) at which
     checkpoint/restart files are created., -1 (no checkpointing)}
   @kw{checkRanges, Enables run-tim checking for membrane voltages that
     are outside of a defined range.  The range is currently hardcoded to
     -110 mV to 60 mV.  A warning will be printed for each cell that has
     voltage outside that range.  Range checking is on by default for
     both the serial and parallel simulation loops.  Set to zero to
     disable the checks., 1}
   @kw{decomposition, The name of the DECOMPOSITION object for this
     simulation., decomposition}
   @kw{diffusion, The name of the DIFFUSION object for this simulation.,
     diffusion}
   @kw{heap, Storage allocated for IO buffers, 500}
   @kw{dt, The time step., 0.01 msec}
   @kw{loop, The initial loop count for the simulation., 0}
   @kw{maxLoop, The maximum value for the loop count., 1000}
   @kw{printRate, , }
   @kw{reaction, The name of the REACTION object for this simulation., reaction}
   @kw{sensor, The name of the sensor object(s) for this simulation.
     Multiple sensors may be specified., No sensors}
   @kw{stateFile, The name of the file(s) from which to load cell model
     state data.  Multiple files may be specified.  They are loaded in
     the order specified.  If the same field is present in multiple
     files\, each load will overwrite previous values., No default.} 
   @kw{stimulus, The name of the STIMULUS object(s) to use in this
     simulation. Multiple stimuli may be specified.  If this keyword is
     not specified there will be no external stimulus in the simulation.}
     {No Stimulus}
   @kw{time, The initial simulation time., 0 msec}
   @endkeywords
*/

void initializeSimulate(const string& name, Simulate& sim)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   sim.name_ = name;

   OBJECT* obj = objectFind(name, "SIMULATE");

   checkForObsoleteKeywords(obj);
   
   int heapSize;
   objectGet(obj, "heap", heapSize, "500");
   heap_start(heapSize);
   
   objectGet(obj, "checkRanges", sim.checkRange_.on, "1");
   objectGet(obj, "VmMin", sim.checkRange_.vMin, "-150");
   objectGet(obj, "VmMax", sim.checkRange_.vMax, " 100");
   objectGet(obj, "loop", (int&)sim.loop_, "0"); // cast away volatile
   objectGet(obj, "maxLoop", sim.maxLoop_, "1000");
   objectGet(obj, "dt", sim.dt_, "0.01", "t");
   objectGet(obj, "time", sim.time_, "0", "t");
   objectGet(obj, "printRate", sim.printRate_, "20");
   objectGet(obj, "globalSyncRate", sim.globalSyncRate_, "-1");
   objectGet(obj, "printGid", sim.printGid_, "-1");
   objectGet(obj, "checkpointRate", sim.checkpointRate_, "-1");
   {
      int tmp; objectGet(obj, "profileAllCounters", tmp, "0");
      if (tmp == 1)
         profileSetVerbosity(true);
      else
         profileSetVerbosity(false);
   }
   {
      int tmp;         objectGet(obj, "writeTorusMap", tmp, "1");
      string filename; objectGet(obj, "torusMapFile", filename, "torusMap");
      if (tmp == 1)
         writeTorusMap(MPI_COMM_WORLD, filename);
   }
   {
      string tmp; objectGet(obj, "checkpointType", tmp, "ascii");
      if (tmp != "ascii")
         sim.asciiCheckpoints_ = false;
      else
         sim.asciiCheckpoints_ = true;
   }
   {
      unsigned nFiles; objectGet(obj, "nFiles", nFiles, "0");
      if (nFiles > 0)
         Pio_setNumWriteFiles(nFiles);
   }
      
   timestampBarrier("initializing anatomy", MPI_COMM_WORLD);
   string nameTmp;
   objectGet(obj, "anatomy", nameTmp, "anatomy");
   initializeAnatomy(sim.anatomy_, nameTmp, MPI_COMM_WORLD);
   sim.nx_ = sim.anatomy_.nx();
   sim.ny_ = sim.anatomy_.ny();
   sim.nz_ = sim.anatomy_.nz();

   objectGet(obj, "stateFile", sim.stateFilename_);
   for (unsigned ii=0; ii<sim.stateFilename_.size(); ++ii)
   {
      string& name = sim.stateFilename_[ii];
      if (name[name.size()-1] != '#')
         name += "#";
   }

   // read either the loopType or parallelDiffusionReaction
   string tmp; objectGet(obj, "loopType", tmp, "omp");
   if (tmp == "omp")
       sim.loopType_ = Simulate::omp;
   else if (tmp == "pdr")
      sim.loopType_ = Simulate::pdr;
   else
      assert(false);
   if (object_testforkeyword(obj, "parallelDiffusionReaction"))
   {
      // backward compatibility
      int tmp; 
      objectGet(obj, "parallelDiffusionReaction", tmp, "0");
      if (tmp == 1)
         sim.loopType_ = Simulate::pdr;
      else
         sim.loopType_ = Simulate::omp;
   }
   timestampBarrier("assigning cells to tasks", MPI_COMM_WORLD);
   string decompositionName;
   objectGet(obj, "decomposition", decompositionName, "decomposition");
   LoadLevel loadLevel = assignCellsToTasks(sim, decompositionName, MPI_COMM_WORLD);

   // default number of diffusion cores is 1 unless the load leveler
   // said otherwise.
   string defaultNDiffusionCores = "1";
   if (loadLevel.nDiffusionCoresHint > 0)
   {
      stringstream stream;
      stream << loadLevel.nDiffusionCoresHint;
      defaultNDiffusionCores  = stream.str();
   }
   
   unsigned nDiffusionCores;
   objectGet(obj, "nDiffusionCores", nDiffusionCores, defaultNDiffusionCores);

   vector<unsigned> diffusionCores;
   objectGet(obj, "diffusionThreads", diffusionCores);

   if (sim.loopType_ == Simulate::pdr)
   {
      // diffusionThreads overrides nDiffusionCores, but when no thread
      // list is specified, we use 1 core.
      if (diffusionCores.size() == 0)
         buildCoreList(nDiffusionCores, diffusionCores);
      ThreadServer& threadServer = ThreadServer::getInstance();
      sim.diffusionThreads_ = threadServer.getThreadTeam(diffusionCores);
      if (getRank(0) == 0)
         cout << "Diffusion Threads: " << sim.diffusionThreads_ << endl;
      sim.reactionThreads_ = threadServer.getThreadTeam(vector<unsigned>());
      if (getRank(0) == 0)
         cout << "Reaction Threads: " << sim.reactionThreads_ << endl;
   }
   
   timestampBarrier("building reaction object", MPI_COMM_WORLD);
   vector<string> reactionNames;
   objectGet(obj, "reaction", reactionNames);
   sim.reaction_ = new ReactionManager;
   for (int ii=0; ii<reactionNames.size(); ++ii) {
      const string& reactionName(reactionNames[ii]);
      sim.reaction_->addReaction(reactionName);
   }
   std::vector<int> cellTypes(sim.anatomy_.nLocal());
   for (int ii=0; ii<cellTypes.size(); ii++)
   {
      cellTypes[ii] = sim.anatomy_.cellType(ii);
   }
   sim.reaction_->create(sim.dt_, cellTypes, sim.reactionThreads_);
   timestampBarrier("finished building reaction object", MPI_COMM_WORLD);

   sim.printIndex_ = -1;
   // -2 -> print index 0 rank 0
   if (sim.printGid_ == -2 && myRank == 0)
   {
      sim.printGid_ = sim.anatomy_.gid(0);
   }
   // -1 -> global min gid
   if ( sim.printGid_ == -1 )    
      sim.printGid_ = findGlobalMinGid(sim.anatomy_);
   
   for (unsigned ii=0; ii<sim.anatomy_.nLocal(); ii++)
      if (sim.anatomy_.gid(ii) == sim.printGid_)
      {
         sim.printIndex_ = ii;
         break;
      }
   
   sim.printFile_=NULL; 
   if (sim.printIndex_ >=0)
   {
      sim.printFile_=fopen("data","a"); 
      printf(                "#   Loop     Time         gid            Vm(t)              dVm_r(t-h)             dVm_d(t-h)\n");
      fprintf(sim.printFile_,"#   Loop     Time         gid            Vm(t)              dVm_r(t-h)             dVm_d(t-h)\n");
   }
   
   
   getRemoteCells(sim, decompositionName, MPI_COMM_WORLD);

   timestampBarrier("building diffusion object", MPI_COMM_WORLD);
   objectGet(obj, "diffusion", nameTmp, "diffusion");
   sim.diffusion_ = diffusionFactory(nameTmp, sim.anatomy_, sim.diffusionThreads_,
                                     sim.reactionThreads_,
                                     sim.loopType_,loadLevel.variantHint);
   
   timestampBarrier("building stimulus object", MPI_COMM_WORLD);
   vector<string> names;
   objectGet(obj, "stimulus", names);
   for (unsigned ii=0; ii<names.size(); ++ii)
   {
      Stimulus* stim = stimulusFactory(names[ii], sim.anatomy_);
      if (stim->nStim() > 0)
	 sim.stimulus_.push_back(stim);
      else
	 delete stim;
   }

   timestampBarrier("building sensor object", MPI_COMM_WORLD);
   names.clear();
   objectGet(obj, "sensor", names);
   for (unsigned ii=0; ii<names.size(); ++ii)
      sim.sensor_.push_back(sensorFactory(names[ii], sim));

}


namespace
{
   void checkForObsoleteKeywords(OBJECT* obj)
   {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
      if (myRank != 0) return;

      if (object_testforkeyword(obj, "snapshotCellList") != 0)
      {
         printf("Obsolete keyword snapshotCellList found in %s object.\n", obj->objclass);
         printf("  Use State Variable Sensor instead\n");
         exit(1);
      }
      if (object_testforkeyword(obj, "snapshotRate") != 0)
      {
         printf("Obsolete keyword snapshotRate found in %s object.\n", obj->objclass);
         printf("  Use State Variable Sensor instead\n");
         exit(1);
      }
      
   }
}

namespace
{
   /** Makes up a core list with one thread per core, except on BGQ
    *  where we have four threads per core. */
   void buildCoreList(unsigned& nCores, vector<unsigned>& cores)
   {
      int threadsPerCore = 1;
      #ifdef BGQ
      threadsPerCore = 4;
      #endif
      assert(cores.size() == 0);
      int nThreads = nCores*threadsPerCore;
      cores.resize(nThreads);
      for (unsigned ii=0; ii<nThreads; ++ii)
         cores[ii] = ii/threadsPerCore;
   }
}

namespace
{
   Long64 findGlobalMinGid(const Anatomy& anatomy)
   {
      Long64 minGid=anatomy.nx()*anatomy.ny()*anatomy.nz();
      for (unsigned ii=0; ii<anatomy.nLocal(); ++ii)
         if (anatomy.gid(ii) < minGid)
            minGid = anatomy.gid(ii); 
      Long64 globalMin;
      MPI_Allreduce(&minGid, &globalMin, 1, MPI_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
      return globalMin;
   }
}


namespace
{
   void writeTorusMap(MPI_Comm comm, const string& filename)
   {
      if ( hi_hasTorus() == 0)
         return;

      int myRank;
      MPI_Comm_rank(comm, &myRank);

      int nDims = hi_nTorusDim();
      int coord[nDims];
      hi_torusCoords(coord);
      
      PFILE* pfile = Popen(filename.c_str(), "w", comm);
      PioSet(pfile, "ngroup", 1);
      if (myRank == 0)
      {
         int size[nDims];
         hi_torusSize(size);
         Pprintf(pfile, "# mpiRank coords[%d]\n", nDims);
         Pprintf(pfile, "# torus size:");
         for (unsigned ii=0; ii<nDims; ++ii)
            Pprintf(pfile, " %4d", size[ii]);
         Pprintf(pfile, "\n");
      }
      
      Pprintf(pfile, "%7d", myRank);
      for (unsigned ii=0; ii<nDims; ++ii)
         Pprintf(pfile, " %4d", coord[ii]);
      Pprintf(pfile, "\n");

      Pclose(pfile);      
   }
}
