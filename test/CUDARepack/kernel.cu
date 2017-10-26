#include <stdio.h>
//#include "types.h"
#include <stdlib.h>

#define CUDA_VERIFY(call) (call)

typedef double Real;

__global__ void setV(long n,Real* tar)
{
  long stride = gridDim.x * blockDim.x;
  long offset = blockDim.x * blockIdx.x + threadIdx.x;
  for(long ii=offset;ii<n;ii+=stride)
  {
     tar[ii] = 0.001 * Real(ii);
  }
}

__global__ void gather(long n,Real* new_data, Real* data, unsigned* index)
{
  int stride = gridDim.x * blockDim.x;
  long offset = blockDim.x * blockIdx.x + threadIdx.x;
  for(long ii=offset;ii<n;ii+=stride)
  {
    new_data[ii] = data[index[ii]];
  }
}

__global__ void scatter(long n,Real* new_data, Real* data, unsigned* index)
{
  int stride = gridDim.x * blockDim.x;
  long offset = blockDim.x * blockIdx.x + threadIdx.x;
  for(long ii=offset;ii<n;ii+=stride)
  {
    new_data[index[ii]] = data[ii];
  }
}

int main(int argc, char** argv)
{
  long num_elem=1000000;

  if (argc == 2)
  {
    num_elem = atoi(argv[1]);
  }

  printf("performing repacking %d elements\n",num_elem);
  srand(123);

  //generating random
  Real *data,*data2,*h_data;
  unsigned *index, *d_index;

  CUDA_VERIFY(cudaMalloc(&data, sizeof(Real)*num_elem));
  CUDA_VERIFY(cudaMalloc(&data2, sizeof(Real)*num_elem));
  CUDA_VERIFY(cudaMalloc(&d_index, sizeof(unsigned)*num_elem));

  h_data = (Real*)malloc(sizeof(Real)*100);
  index = (unsigned*)malloc(sizeof(unsigned)*num_elem);

  for(int ii=0;ii<num_elem;ii++) index[ii]=ii;

  for(int ii=0;ii<num_elem*10;ii++)
  {
    unsigned r1 = rand();
    unsigned r2 = rand();

    unsigned rr = r1 * r2;

    rr = rr % num_elem;
  
    unsigned tmp = index[rr];
    index[rr] = index[ii % num_elem];
    index[ii % num_elem] = tmp;
  }

  CUDA_VERIFY(cudaMemcpy(d_index,index,sizeof(unsigned)*num_elem,cudaMemcpyHostToDevice));

  setV<<<112,1024>>>(num_elem,data);

  gather<<<224,1024>>>(num_elem,data2,data,d_index);
  scatter<<<224,1024>>>(num_elem,data2,data,d_index);

  CUDA_VERIFY(cudaMemcpy(h_data,data2,sizeof(Real)*100,cudaMemcpyDeviceToHost));

  for(int ii=0;ii<100;ii++) printf("%d=%f\n",ii,h_data[ii]);

  free(h_data);

  
  CUDA_VERIFY(cudaFree(data));
  CUDA_VERIFY(cudaFree(data2));
  CUDA_VERIFY(cudaFree(d_index));

}



