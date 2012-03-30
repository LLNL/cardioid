#include "initializeSimulate.hh"
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
#include "Threading.hh"
#include "mpiUtils.h"
#include "PerformanceTimers.hh"

using namespace std;

void initializeSimulate(const string& name, Simulate& sim)
{
   sim.name_ = name;

   OBJECT* obj = objectFind(name, "SIMULATE");

   objectGet(obj, "loop", sim.loop_, "0");
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
   objectGet(obj, "nDiffusionCores", nDiffusionCores, "1");
   // Need the diffusionGroup to be initialized for loopIO in simulationLoop.cc
   sim.diffusionGroup_ = sim.tinfo_.mkGroup(nDiffusionCores); 
   if (sim.parallelDiffusionReaction_ == 1)    
   {
      sim.reactionGroup_  = sim.tinfo_.mkGroup(-1); 
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
   objectGet(obj, "decomposition", nameTmp, "decomposition");
   assignCellsToTasks(sim, nameTmp, MPI_COMM_WORLD);
   
   getRemoteCells(sim, nameTmp, MPI_COMM_WORLD);

   timestampBarrier("building diffusion object", MPI_COMM_WORLD);
   objectGet(obj, "diffusion", nameTmp, "diffusion");
   sim.diffusion_ = diffusionFactory(nameTmp, sim.anatomy_, *sim.diffusionGroup_,
                                     sim.parallelDiffusionReaction_);
   
   timestampBarrier("building reaction object", MPI_COMM_WORLD);
   objectGet(obj, "reaction", nameTmp, "reaction");
   sim.reaction_ = reactionFactory(nameTmp, sim.anatomy_,sim.reactionGroup_);

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

}
