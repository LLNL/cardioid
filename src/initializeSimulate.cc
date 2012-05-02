#include "initializeSimulate.hh"

#include <iostream>

#include "object_cc.hh"
#include "Simulate.hh"

#include "initializeAnatomy.hh"
#include "assignCellsToTasks.hh"
#include "diffusionFactory.hh"
#include "reactionFactory.hh"
#include "stimulusFactory.hh"
#include "sensorFactory.hh"
#include "getRemoteCells.hh"
#include "Anatomy.hh"
#include "mpiUtils.h"
#include "PerformanceTimers.hh"
#include "readSnapshotCellList.hh"

using namespace std;


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


void initializeSimulate(const string& name, Simulate& sim)
{
   sim.name_ = name;

   OBJECT* obj = objectFind(name, "SIMULATE");

   {
     /*
       Since sim.loop_ is volatile, objectGet(...)
       does not accept it as an argument of type int &;
       hence the workaround with a temporary variable.
     */
     int loopInit;
     objectGet(obj, "loop", loopInit, "0");
     sim.loop_ = loopInit;
   }
   objectGet(obj, "maxLoop", sim.maxLoop_, "1000");
   objectGet(obj, "dt", sim.dt_, "0.01", "t");
   objectGet(obj, "time", sim.time_, "0", "t");
   objectGet(obj, "printRate", sim.printRate_, "20");
   objectGet(obj, "snapshotRate", sim.snapshotRate_, "100");
   objectGet(obj, "checkpointRate", sim.checkpointRate_, "-1");
   objectGet(obj, "parallelDiffusionReaction", sim.parallelDiffusionReaction_, "0");
   unsigned profileAllCounters_;
   objectGet(obj, "profileAllCounters", profileAllCounters_, "0");
   if (profileAllCounters_ == 1)
      profileSetVerbosity(true);
   else
      profileSetVerbosity(false);
   unsigned nDiffusionCores;
   vector<unsigned> diffusionCores;
   objectGet(obj, "nDiffusionCores", nDiffusionCores, "1");
   objectGet(obj, "diffusionThreads", diffusionCores);

   if (sim.parallelDiffusionReaction_ == 1)
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
   timestampBarrier("assigning cells to tasks", MPI_COMM_WORLD);
   string decompositionName;
   objectGet(obj, "decomposition", decompositionName, "decomposition");
   assignCellsToTasks(sim, decompositionName, MPI_COMM_WORLD);
   
   timestampBarrier("building reaction object", MPI_COMM_WORLD);
   objectGet(obj, "reaction", nameTmp, "reaction");
   sim.reaction_ = reactionFactory(nameTmp, sim.anatomy_, sim.reactionThreads_);

   timestampBarrier("getRemoteCells", MPI_COMM_WORLD);
   getRemoteCells(sim, decompositionName, MPI_COMM_WORLD);

   timestampBarrier("building diffusion object", MPI_COMM_WORLD);
   objectGet(obj, "diffusion", nameTmp, "diffusion");
   sim.diffusion_ = diffusionFactory(nameTmp, sim.anatomy_, sim.diffusionThreads_,
                                     sim.reactionThreads_,
                                     sim.parallelDiffusionReaction_);
   
   timestampBarrier("building stimulus object", MPI_COMM_WORLD);
   vector<string> names;
   objectGet(obj, "stimulus", names);
   for (unsigned ii=0; ii<names.size(); ++ii)
      sim.stimulus_.push_back(stimulusFactory(names[ii], sim.anatomy_));

   timestampBarrier("building sensor object", MPI_COMM_WORLD);
   names.clear();
   objectGet(obj, "sensor", names);
   for (unsigned ii=0; ii<names.size(); ++ii)
      sim.sensor_.push_back(sensorFactory(names[ii], sim.anatomy_));

   // let user specify a filename containing the list of cell gids they want in the snapshots
   // (need to do this after data is distributed, so we can store just the local subset of points
   // on each task)
   string snapshotCLFile;
   objectGet(obj, "snapshotCellList", snapshotCLFile, "");
   if (snapshotCLFile != "")
      readSnapshotCellList(snapshotCLFile, sim,obj);
}
