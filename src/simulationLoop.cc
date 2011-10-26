#include "simulationLoop.hh"

#include <vector>

#include "Simulate.hh"
#include "Diffusion.hh"
#include "Reaction.hh"
#include "HaloExchange.hh"
#include "GridRouter.hh"

using namespace std;

void simulationLoop(Simulate& sim)
{
  
   vector<double> iStimArray(sim.anatomy_.nLocal());

   sim.voltageExchange_ = new HaloExchange<double>(sim.router_->sendMap(), sim.router_->commTable());
   
   for ( ; sim.loop_<sim.maxLoop_; ++sim.loop_)
   {
      sim.voltageExchange_->execute(sim.VmArray_, sim.anatomy_.nLocal());
      
      // DIFFUSION
      sim.diffusion_->diffusion(sim.VmArray_, iStimArray);
      
      // REACTION
// 	 iStimArray[ii] *= param.diffusionscale;
      
// 	 // code to limit or set iStimArray goes here.
      
      sim.reaction_->calc(sim.dt_, sim.VmArray_, iStimArray);
   }
}
