#include "Ledger.hh"

__global__ void fillSendBuff_kernel(double* sendBufRaw, double* dataRaw, int* sendMapRaw, const int begin, const int end)
{
   int ii = threadIdx.x + blockIdx.x*blockDim.x  +begin;
   if (ii >= end) { return; }
   sendBufRaw[ii] = dataRaw[sendMapRaw[ii]];
}


void fillSendBufferCUDA(const double* sendBufRaw, const double* dataRaw, const int* sendMapRaw, const int begin, const int end){
   int blockSize = 1024;
   fillSendBuff_kernel<<<(end-begin+blockSize-1)/blockSize,blockSize>>>
      (ledger_lookup(sendBufRaw),
       ledger_lookup(dataRaw),
       ledger_lookup(sendMapRaw),
       begin,
       end);
}
