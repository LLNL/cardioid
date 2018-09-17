#include "lazy_array.hh"

__global__ void fillSendBuff_kernel(wo_larray_ptr<double> sendBufRaw, 
                        ro_larray_ptr<double> dataRaw, 
                        ro_larray_ptr<int> sendMapRaw,
                        const int begin, const int end)
{
   int ii = threadIdx.x + blockIdx.x*blockDim.x  +begin;
   if (ii >= end) { return; }
   sendBufRaw[ii] = dataRaw[sendMapRaw[ii]];
}


void fillSendBufferCUDA(wo_larray_ptr<double> sendBufRaw, 
                        ro_larray_ptr<double> dataRaw, 
                        ro_larray_ptr<int> sendMapRaw)
{
   ContextRegion region(GPU);
   int blockSize = 1024;
   int begin=0;
   int end=sendMapRaw.size();
   if (sendMapRaw.size() > 0)
   {
      fillSendBuff_kernel<<<(end-begin+blockSize-1)/blockSize,blockSize>>>
         (sendBufRaw,
          dataRaw,
          sendMapRaw,
          begin,
          end);
   }
}
