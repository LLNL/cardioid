
#include "lazy_array.hh"

__global__ void addStimulusKernel(rw_larray_ptr<double> dVmDiffusion, ro_larray_ptr<int> stimList, const double value)
{
   int ii = threadIdx.x + blockIdx.x*blockDim.x;
   if (ii >= stimList.size()) { return; }
   dVmDiffusion[stimList[ii]] += value;
}

void addStimulus(rw_larray_ptr<double> dVmDiffusion, ro_larray_ptr<int> stimList, const double value)
{
   ContextRegion region(GPU);
   int blockSize=1024;
   addStimulusKernel<<<(stimList.size()+blockSize-1)/blockSize,blockSize>>>(
      dVmDiffusion,
      stimList,
      value);
}
      
