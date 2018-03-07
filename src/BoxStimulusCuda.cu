
#include "TransportCoordinator.hh"

__global__ void addStimulusKernel(ArrayView<double> dVmDiffusionRaw, ConstArrayView<int> stimListRaw, const double value)
{
   int ii = threadIdx.x + blockIdx.x*blockDim.x;
   if (ii >= stimListRaw.size()) { return; }
   dVmDiffusionRaw[stimListRaw[ii]] += value;
}

void addStimulus(OnDevice<ArrayView<double>> dVmDiffusion, OnDevice<ConstArrayView<int>> stimList, const double value)
{
   int blockSize=1024;
   addStimulusKernel<<<(stimList.size()+blockSize-1)/blockSize,blockSize>>>(
      dVmDiffusion,
      stimList,
      value);
}
      
