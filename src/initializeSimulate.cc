#include <omp.h> 
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
   sim.diffusionGroup_ = sim.tinfo_.mkGroup(1); 
   sim.reactionGroup_  = sim.tinfo_.mkGroup(-1); 
   

   string nameTmp;
   objectGet(obj, "anatomy", nameTmp, "anatomy");
   initializeAnatomy(sim.anatomy_, nameTmp, MPI_COMM_WORLD);
   sim.nx_ = sim.anatomy_.nx();
   sim.ny_ = sim.anatomy_.ny();
   sim.nz_ = sim.anatomy_.nz();

   objectGet(obj, "stateFile", sim.stateFilename_, "");
   if (!sim.stateFilename_.empty())
      if (sim.stateFilename_[sim.stateFilename_.size()-1] != '#')
          sim.stateFilename_ += "#";
   objectGet(obj, "decomposition", nameTmp, "decomposition");
   assignCellsToTasks(sim, nameTmp, MPI_COMM_WORLD);
   
   getRemoteCells(sim, nameTmp, MPI_COMM_WORLD);

// Assume thread are assigned in round robin order. 

   objectGet(obj, "diffusion", nameTmp, "diffusion");
   sim.diffusion_ = diffusionFactory(nameTmp, sim.anatomy_);
   
   objectGet(obj, "reaction", nameTmp, "reaction");
   sim.reaction_ = reactionFactory(nameTmp, sim.anatomy_,sim.reactionGroup_);

   vector<string> names;
   objectGet(obj, "stimulus", names);
   for (unsigned ii=0; ii<names.size(); ++ii) sim.stimulus_.push_back(stimulusFactory(names[ii], sim.anatomy_));

   names.clear();
   objectGet(obj, "sensor", names);
   for (unsigned ii=0; ii<names.size(); ++ii) sim.sensor_.push_back(sensorFactory(names[ii], sim.anatomy_));

}
