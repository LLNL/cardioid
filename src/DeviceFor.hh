#pragma once

#include "lazy_array.hh"

#define HOST_PARALLEL_FORALL(N, VAR, BODY)      \
do                                              \
{                                               \
   for (int VAR=0; VAR<(N); VAR++)              \
   {                                            \
      BODY;                                     \
   }                                            \
} while(0)

#ifdef USE_CUDA

template<typename LoopBody>
__global__ void forall_device_kernel(const int len, LoopBody loopBody)
{
   const int _ii = threadIdx.x + blockIdx.x*blockDim.x;
   if (_ii >= len) { return; }
   loopBody(_ii);
}

#define DEVICE_PARALLEL_FORALL(N, VAR, LAMBDA) \
do                                                                      \
{                                                                       \
   int blockSize = 128;                                                 \
   forall_device_kernel<<<(N+blockSize-1)/blockSize,blockSize>>>(N, [=] __device__ (const int VAR) mutable {LAMBDA;}); \
} while (0)


#else

#define DEVICE_PARALLEL_FORALL(N, VAR, LAMBDA) HOST_PARALLEL_FORALL(N, VAR, LAMBDA)

#endif
