#pragma once

#include "lazy_array.hh"

#define HOST_PARALLEL_FORALL(N, LAMBDA)         \
do                                              \
{                                               \
   ContextRegion region(CPU);                   \
   auto lambda = [=] LAMBDA;                    \
   for (int ii=0; ii<(N); ii++)                 \
   {                                            \
      lambda(ii);                               \
   }                                            \
} while(0)

#ifdef USE_CUDA

#else

#define DEVICE_PARALLEL_FORALL(N, LAMBDA) HOST_PARALLEL_FORALL(N, LAMBDA)

#endif
