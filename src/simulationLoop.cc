#include "simulationLoop.hh"

#include <vector>

#include "Simulate.hh"
#include "Diffusion.hh"
#include "Reaction.hh"

using namespace std;

void simulationLoop(Simulate& sim)
{
  
   vector<double> iStimArray(sim.nCellLocal_);
   
   for ( ; sim.loop_<sim.maxLoop_; ++sim.loop_)
   {
      
      // DIFFUSION
      sim.diffusion_->diffusion(sim.VmArray_, iStimArray);
      
      // REACTION
// 	 iStimArray[ii] *= param.diffusionscale;
      
// 	 // code to limit or set iStimArray goes here.
      
      sim.reaction_->calc(sim.dt_, sim.VmArray_, iStimArray);
   }
}
