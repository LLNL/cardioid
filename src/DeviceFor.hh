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

#define DEVICE_PARALLEL_FORALL(N, VAR, LAMBDA) HOST_PARALLEL_FORALL(N, VAR, LAMBDA)

#else

#define DEVICE_PARALLEL_FORALL(N, VAR, LAMBDA) HOST_PARALLEL_FORALL(N, VAR, LAMBDA)

#endif
