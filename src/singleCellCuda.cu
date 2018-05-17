
#include "TransportCoordinator.hh"

__global__ void setStimulusKernel(ArrayView<double> iStim, const double stimAmount) {
   int ii = threadIdx.x + blockIdx.x*blockDim.x;
   if (ii >= iStim.size()) { return; }
   //use a negative sign here to undo the negative we had above.
   iStim[ii] = -stimAmount;
}

void setStimulus(OnDevice<ArrayView<double>> iStim, const double stimAmount) {
   int blockSize = 1024;
   int nCells = iStim.size();
   
   setStimulusKernel<<<(nCells+blockSize-1)/blockSize,blockSize>>>(
      iStim,
      stimAmount);
}
