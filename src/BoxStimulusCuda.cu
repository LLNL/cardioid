
#include "Ledger.hh"

__global__ void addStimulusKernel(double* dVmDiffusionRaw, const int* stimListRaw, const int stimListSize, const double value)
{
   int ii = threadIdx.x + blockIdx.x*blockDim.x;
   if (ii >= stimListSize) { return; }
   dVmDiffusionRaw[stimListRaw[ii]] += value;
}

void addStimulus(double* dVmDiffusionRaw, const int* stimListRaw, const int stimListSize, const double value)
{
   int blockSize=1024;
   addStimulusKernel<<<(stimListSize+blockSize-1)/blockSize,blockSize>>>(
      ledger_lookup(dVmDiffusionRaw),
      ledger_lookup(stimListRaw),
      stimListSize,
      value);
}
      
