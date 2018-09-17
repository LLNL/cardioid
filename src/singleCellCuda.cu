
#include "lazy_array.hh"

__global__ void setStimulusKernel(wo_larray_ptr<double> iStim, const double stimAmount) {
   int ii = threadIdx.x + blockIdx.x*blockDim.x;
   if (ii >= iStim.size()) { return; }
   //use a negative sign here to undo the negative we had above.
   iStim[ii] = -stimAmount;
}

void setStimulus(wo_larray_ptr<double> iStim, const double stimAmount) {
   int blockSize = 1024;
   int nCells = iStim.size();
   ContextRegion region(GPU);

   setStimulusKernel<<<(nCells+blockSize-1)/blockSize,blockSize>>>(
      iStim,
      stimAmount);
}
