#include "simulationLoopCuda.hh"

__global__ void setStimulusKernel(ArrayView<double> iStimRaw,ConstArrayView<double> dVmDiffusionRaw)
{
   int ii = threadIdx.x + blockIdx.x*blockDim.x;
   if (ii >= iStimRaw.size()) { return; }
   iStimRaw[ii] = -dVmDiffusionRaw[ii];
}

void setStimulus(OnDevice<ArrayView<double>> iStim,OnDevice<ConstArrayView<double>> dVmDiffusionRaw)
{
   int blockSize = 1024;
   setStimulusKernel<<<(iStim.size()+blockSize-1)/blockSize,blockSize>>>
      (iStim,
       dVmDiffusionRaw);
}


__global__ void integrateVoltageKernel(ArrayView<double> vmarrayRaw, ConstArrayView<double> dVmReactionRaw, ConstArrayView<double> dVmDiffusionRaw, const double dt)
{
   int ii = threadIdx.x + blockIdx.x*blockDim.x;
   if (ii >= dVmDiffusionRaw.size()) { return; }
   vmarrayRaw[ii] += dt*(dVmReactionRaw[ii] + dVmDiffusionRaw[ii]);
}

void integrateVoltage(OnDevice<ArrayView<double>> vmarrayRaw, OnDevice<ConstArrayView<double>> dVmReactionRaw, OnDevice<ConstArrayView<double>> dVmDiffusionRaw,const double dt)
{
   int blockSize = 1024;
   integrateVoltageKernel<<<(dVmDiffusionRaw.size()+blockSize-1)/blockSize,blockSize>>>
      (vmarrayRaw,
       dVmReactionRaw,
       dVmDiffusionRaw,
       dt);
}
