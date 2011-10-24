#include "simulationLoop.hh"

#include <vector>

#include "Simulate.hh"
#include "Diffusion.hh"

using namespace std;

/** The BlueBeats code uses the following data structures:
 *
 *  LINEAR ARRAYS
 *
 *  VmArray     voltage
 *  IstimArray
 *  pemIBMArray The cell models
 *  cells       Array of tuples to work on
 *
 *  3D ARRAYS
 *
 *  Vm
 *  diffIntra
 *
 *
 *  MAPPING ARRAYS
 *
 *  VM3DPointer     Used at the end of the time step to copy the 
 *                  values from VmArray into Vm
*/


/**
   

*/
void simulationLoop(Simulate& sim)
{
  
   vector<double> iStimArray(sim.nCellLocal_);
   
   for ( ; sim.loop_<sim.maxLoop_; ++sim.loop_)
   {
      
      // DIFFUSION
      sim.diffusion_->diffusion(sim.VmArray_, iStimArray);
      
//       // REACTION
//       for (int iCell=0; iCell<nTissue; ++iCell)
//       {
// 	 iStimArray[iCell] *= param.diffusionscale;
      
// 	 // code to limit or set iStimArray goes here.
      
// 	 VmArray[iCell] = pemIBMArray[iCell]->Calc(
// 	    param.dt, VmArray[iCell], IstimArray[iCell]);
//       }
   }
}
