#include "Ledger.hh"
#include "TransportCoordinator.hh"
#include "Long64.hh"
#include <stdio.h>

#define NUM_TBLOCK 112
#define CELL_PER_BLOCK 32


__global__ void calcInvr_kernel(ArrayView<double> invr,
                                ConstArrayView<Long64> gids,
                                ConstArrayView<double> ecgPoints,
                                const int nEcgPoints,
                                const int nx, const int ny, const int nz,
                                const double dx, const double dy, const double dz,
                                const int begin, const int end)
{

   int ii = threadIdx.x + blockIdx.x*blockDim.x  +begin;
   if (ii >= end) { return; }

    const int dim=3;

    Long64 gid=gids[ii];
    int x=gid%nx;
    int y=(gid/nx) %ny;
    int z=gid/nx/ny;

    double xcoor=x*dx;
    double ycoor=y*dy;
    double zcoor=z*dz;

    for(int j=0; j<nEcgPoints; ++j)
    {
        double dxx=xcoor-ecgPoints[j*dim];
        double dyy=ycoor-ecgPoints[j*dim+1];
        double dzz=zcoor-ecgPoints[j*dim+2];
        invr[ii*nEcgPoints+j]=1.0/sqrt(dxx*dxx+dyy*dyy+dzz*dzz);
    }

}

void calcInvrCUDA(OnDevice<ArrayView<double>> invr,
                  OnDevice<ConstArrayView<Long64>> gids, 
                  OnDevice<ConstArrayView<double>> ecgPoints,
                  const int nEcgPoints,
                  const int nx, const int ny, const int nz,
                  const double dx, const double dy, const double dz)
{
    int blockSize = 1024;
    ConstArrayView<Long64> tmp=gids;

    int begin=0;
    int end=tmp.size();


    calcInvr_kernel<<<(end-begin+blockSize-1)/blockSize, blockSize>>>
        (invr,
         gids,
         ecgPoints,
         nEcgPoints,
         nx, ny, nz,
         dx, dy, dz,
         begin, end);

}

#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 600
inline __device__ double atomicAdd(double* pointer, const double value)
{
   unsigned long long int* address_as_ull = reinterpret_cast<unsigned long long int*>(pointer);
   unsigned long long int old = *address_as_ull;
   unsigned long long int assumed;
   do {
      assumed=old;
      old = atomicCAS(address_as_ull,
                      assumed,
                      __double_as_longlong(value + __longlong_as_double(assumed)));
   } while(assumed != old);
   return assumed;
}
#endif

__global__ void calcEcg_kernel(ArrayView<double> ecgs,
                               ConstArrayView<double> invr,
                               ConstArrayView<double> dVmDiffusion,
                               const int nEcgPoints,
                               const int cellPartition,
                               const int nCells)
{
  __shared__ double smResult[1024];
  double tmpECG=0;

  const int cellStart = threadIdx.x / nEcgPoints + blockIdx.x * cellPartition;
  const int cellStride = blockDim.x / nEcgPoints;
  const int cellEnd = (blockIdx.x<(NUM_TBLOCK-1))?((blockIdx.x+1)*cellPartition):(nCells);
  const int ecgID = threadIdx.x % nEcgPoints;

  const unsigned int ecgSet = blockDim.x/ nEcgPoints;    

  if ( threadIdx.x >= ecgSet * nEcgPoints ) return;

  for(int cell=cellStart;cell<cellEnd;cell += cellStride)
  {
    tmpECG +=dVmDiffusion[cell] * invr[ecgID + nEcgPoints*cell];
  }

  const unsigned int close2N = 0x80000000 >> __clz(ecgSet-1);

  if ( threadIdx.x >= close2N*nEcgPoints )
      smResult[threadIdx.x] = tmpECG;
  __syncthreads();
  if (threadIdx.x < (ecgSet - close2N)*nEcgPoints )
      tmpECG += smResult[threadIdx.x + close2N*nEcgPoints];

  int jump = close2N >> 1;
  while(jump>0)
  {
    if (( threadIdx.x >= jump*nEcgPoints ) && ( threadIdx.x < 2*jump*nEcgPoints ))
      smResult[threadIdx.x] = tmpECG;
    __syncthreads();
    if(threadIdx.x < jump*nEcgPoints)
      tmpECG += smResult[threadIdx.x + jump*nEcgPoints];
    jump = jump >> 1;
  }

  if(threadIdx.x<nEcgPoints){
    double* result=&ecgs[0];
    atomicAdd(result+ecgID ,tmpECG);
   }
   
}

void calcEcgCUDA(OnDevice<ArrayView<double>> ecgs,
                 OnDevice<ConstArrayView<double>> invr,
                 OnDevice<ConstArrayView<double>> dVmDiffusion,
                 const int nEcgPoints)
{

    ConstArrayView<double> tmp=dVmDiffusion;
    int nCells=tmp.size();
    const int cellPartition = (nCells+(NUM_TBLOCK-1))/(NUM_TBLOCK);

    if ( nEcgPoints > 32 ) printf("error:too many ECG points\n");
    calcEcg_kernel<<<NUM_TBLOCK, (nEcgPoints*CELL_PER_BLOCK)>>>
        (ecgs,
         invr,
         dVmDiffusion,
         nEcgPoints,
         cellPartition,
         nCells);
}

__global__ void dump_kernel(ConstArrayView<double> data)
{
	for(int ii=0; ii<100; ii++)
        {
	    printf("%f ",data[ii]);
	    if((ii+1)%10==0) printf("\n");
        }
}

void dump_GPU_data(OnDevice<ArrayView<double>> ecgs,
                 OnDevice<ConstArrayView<double>> invr,
                 OnDevice<ConstArrayView<double>> dVmDiffusion,
                 const int nEcgPoints)
{

    dump_kernel<<<1,1>>>(invr);
    printf("\n\n\n");
    dump_kernel<<<1,1>>>(dVmDiffusion);
    printf("\n\n\n");
}
