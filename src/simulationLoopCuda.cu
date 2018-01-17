#include "simulationLoopCuda.hh"
#include "Ledger.hh"

__global__ void setStimulusKernel(double* iStimRaw,const double* dVmDiffusionRaw, const int nLocal)
{
   int ii = threadIdx.x + blockIdx.x*blockDim.x;
   if (ii >= nLocal) { return; }
   iStimRaw[ii] = -dVmDiffusionRaw[ii];
}

void setStimulus(double* iStimRaw,const double* dVmDiffusionRaw, const int nLocal)
{
   int blockSize = 1024;
   setStimulusKernel<<<(nLocal+blockSize-1)/blockSize,blockSize>>>
      ((double*)ledger_lookup(iStimRaw),
       (double*)ledger_lookup(dVmDiffusionRaw),
       nLocal);
}


__global__ void integrateVoltageKernel(double* vmarrayRaw, const double* dVmReactionRaw, const double* dVmDiffusionRaw, const int nLocal, const double dt)
{
   int ii = threadIdx.x + blockIdx.x*blockDim.x;
   if (ii >= nLocal) { return; }
   vmarrayRaw[ii] += dt*(dVmReactionRaw[ii] + dVmDiffusionRaw[ii]);
}

void integrateVoltage(double* vmarrayRaw, const double* dVmReactionRaw, const double* dVmDiffusionRaw, const int nLocal, const double dt)
{
   int blockSize = 1024;
   integrateVoltageKernel<<<(nLocal+blockSize-1)/blockSize,blockSize>>>
      ((double*)ledger_lookup(vmarrayRaw),
       (double*)ledger_lookup(dVmReactionRaw),
       (double*)ledger_lookup(dVmDiffusionRaw),
       nLocal, dt);
}

