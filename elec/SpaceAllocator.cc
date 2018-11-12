
#include <cstring>
#include <iostream>
#include <cassert>

#ifdef USE_CUDA
#include <cuda.h>
#include <cuda_runtime_api.h>
#define CUDA_VERIFY(x) do { cudaError_t error = x; if (error != cudaSuccess) { std::cout << error << std::endl; assert(error == cudaSuccess && #x ); } } while(0)
#else
#include <cstdlib>
#endif


#include "SpaceAllocator.hh"


template <>
bool spaceMalloc<void>(const ExecutionSpace space,void** dst, const std::size_t size)
{
#ifdef USE_CUDA
   if (space == CPU)
   {
      CUDA_VERIFY(cudaMallocHost(dst, size));
   }
   else
   {
      CUDA_VERIFY(cudaMalloc(dst, size));
   }
   return true;
#else
   return posix_memalign(dst,512,size) != 0;
#endif
}

template <>
void spaceMemcpy<void>(const ExecutionSpace dstSpace, void* dst,
                       const ExecutionSpace srcSpace, const void* src,
                       const std::size_t size)
{
   if (0) {}
#ifdef USE_CUDA
   else if (dstSpace == GPU && srcSpace == GPU)
   {
      CUDA_VERIFY(cudaMemcpy(dst, src, size, cudaMemcpyDeviceToDevice));
   }
   else if (dstSpace == CPU && srcSpace == GPU)
   {
      CUDA_VERIFY(cudaMemcpy(dst, src, size, cudaMemcpyDeviceToHost));
   }
   else if (dstSpace == GPU && srcSpace == CPU)
   {
      CUDA_VERIFY(cudaMemcpy(dst, src, size, cudaMemcpyHostToDevice));
   }
#endif
   else
   {
      memcpy(dst,src,size);
   }
}

template <>
void spaceFree<void>(const ExecutionSpace space,void* dst)
{
#ifdef USE_CUDA
   if (space == CPU)
   {
      CUDA_VERIFY(cudaFreeHost(dst));
   }
   else
   {
      CUDA_VERIFY(cudaFree(dst));
   }
#else
   free(dst);
#endif
}
                   
