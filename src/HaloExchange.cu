#include "Ledger.hh"
#include "TransportCoordinator.hh"

__global__ void fillSendBuff_kernel(ArrayView<double> sendBufRaw, ConstArrayView<double> dataRaw, ConstArrayView<int> sendMapRaw, const int begin, const int end)
{
   int ii = threadIdx.x + blockIdx.x*blockDim.x  +begin;
   if (ii >= end) { return; }
   sendBufRaw[ii] = dataRaw[sendMapRaw[ii]];
}


void fillSendBufferCUDA(OnDevice<ArrayView<double>> sendBufRaw, 
                        OnDevice<ConstArrayView<double>> dataRaw, 
                        OnDevice<ConstArrayView<int>> sendMapRaw)
{
   int blockSize = 1024;
   ConstArrayView<int> temp = sendMapRaw;
   int begin=0;
   int end=temp.size();
   fillSendBuff_kernel<<<(end-begin+blockSize-1)/blockSize,blockSize>>>
      (sendBufRaw,
       dataRaw,
       sendMapRaw,
       begin,
       end);
}
