#include "simulationLoopCuda.hh"

__global__ void setStimulusKernel(wo_larray_ptr<double> iStim,ro_larray_ptr<double> dVmDiffusion)
{
   int ii = threadIdx.x + blockIdx.x*blockDim.x;
   if (ii >= iStim.size()) { return; }
   iStim[ii] = -dVmDiffusion[ii];
}

void setStimulus(wo_larray_ptr<double> iStim,ro_larray_ptr<double> dVmDiffusion)
{
   ContextRegion region(GPU);
   int blockSize = 1024;
   setStimulusKernel<<<(iStim.size()+blockSize-1)/blockSize,blockSize>>>
      (iStim,
       dVmDiffusion);
}


__global__ void integrateVoltageKernel(rw_larray_ptr<double> vmarray,
                                       ro_larray_ptr<double> dVmReaction,
                                       ro_larray_ptr<double> dVmDiffusion,
                                       const double dt)
{
   int ii = threadIdx.x + blockIdx.x*blockDim.x;
   if (ii >= dVmDiffusion.size()) { return; }
   vmarray[ii] += dt*(dVmReaction[ii] + dVmDiffusion[ii]);
}

void integrateVoltage(rw_larray_ptr<double> vmarray,
                      ro_larray_ptr<double> dVmReaction,
                      ro_larray_ptr<double> dVmDiffusion,
                      const double dt)
{
   ContextRegion region(GPU);
   int blockSize = 1024;
   integrateVoltageKernel<<<(dVmReaction.size()+blockSize-1)/blockSize,blockSize>>>
      (vmarray,
       dVmReaction,
       dVmDiffusion,
       dt);
}
