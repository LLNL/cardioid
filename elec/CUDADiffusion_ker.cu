//#include "CUDADiffusion.hh"
//#include "DiffusionUtils.hh"
//#include "SymmetricTensor.hh"
//#include <vector>
//#include <map>

//#include "options.h"
//#include "cudautil.h"
#include <stdio.h>
#include "lazy_array.hh"

#define XTILE 20
typedef double Real;


__global__ void diff_6face_v1(const Real* d_psi, Real* d_npsi, const Real* d_sigmaX, const Real* d_sigmaY, const Real* d_sigmaZ,int Lii, int Ljj, int Lkk)
{

  //map z dir to threads
  //z is the fastest varying direction

  //2d decomposition
  //32x32 in y z direction
  __shared__ Real sm_psi[4][32][32]; //32 KB

  #define V0(y,z) sm_psi[pii][y][z]
  #define V1(y,z) sm_psi[cii][y][z]
  #define V2(y,z) sm_psi[nii][y][z]

  #define sigmaX(x,y,z,dir) d_sigmaX[ z + Lkk * ( y + Ljj * ( x + Lii * dir ) ) ]
  #define sigmaY(x,y,z,dir) d_sigmaY[ z + Lkk * ( y + Ljj * ( x + Lii * dir ) ) ]
  #define sigmaZ(x,y,z,dir) d_sigmaZ[ z + Lkk * ( y + Ljj * ( x + Lii * dir ) ) ]

  #define psi(x,y,z) d_psi[ z + Lkk * ( (y) + Ljj * (x) ) ]
  #define npsi(x,y,z) d_npsi[ z + Lkk * ( (y) + Ljj * (x) ) ]

  int tjj = threadIdx.y;
  int tkk = threadIdx.x;

  //shift for each tile
  d_psi    = &(psi(XTILE*blockIdx.x, 30*blockIdx.y, 30*blockIdx.z));
  d_npsi   = &(npsi(XTILE*blockIdx.x, 30*blockIdx.y, 30*blockIdx.z));

  d_sigmaX  = &(sigmaX(XTILE*blockIdx.x, 30*blockIdx.y, 30*blockIdx.z, 0));
  d_sigmaY  = &(sigmaY(XTILE*blockIdx.x, 30*blockIdx.y, 30*blockIdx.z, 0));
  d_sigmaZ  = &(sigmaZ(XTILE*blockIdx.x, 30*blockIdx.y, 30*blockIdx.z, 0));

  int Last_x=XTILE+1; int nLast_y=31; int nLast_z=31;
  if (blockIdx.x == gridDim.x-1) Last_x = Lii-2 - XTILE * blockIdx.x + 1;
  if (blockIdx.y == gridDim.y-1) nLast_y = Ljj-2 - 30 * blockIdx.y + 1;
  if (blockIdx.z == gridDim.z-1) nLast_z = Lkk-2 - 30 * blockIdx.z + 1;

  if(tjj>nLast_y) return;
  if(tkk>nLast_z) return;

  int pii,cii,nii,tii;
  pii=0; cii=1; nii=2;

  sm_psi[cii][tjj][tkk] = psi(0,tjj,tkk);
  sm_psi[nii][tjj][tkk] = psi(1,tjj,tkk);
  Real xcharge,ycharge,zcharge,dV;

  __syncthreads();
  //initial
  if ((tkk>0) && (tkk<nLast_z) && (tjj>0) && (tjj<nLast_y))
  {
    Real xd=-V1(tjj,tkk) + V2(tjj,tkk);
    Real yd=(-V1(-1 + tjj,tkk) + V1(1 + tjj,tkk) - V2(-1 + tjj,tkk) + V2(1 + tjj,tkk))/4.;
    Real zd=(-V1(tjj,-1 + tkk) + V1(tjj,1 + tkk) - V2(tjj,-1 + tkk) + V2(tjj,1 + tkk))/4.;

    dV = 0;
    dV -= sigmaX(0+1,tjj,tkk,0) * xd + sigmaX(0+1,tjj,tkk,1) * yd + sigmaX(0+1,tjj,tkk,2) * zd ; 
  }

  tii=pii; pii=cii; cii=nii; nii=tii;

  for(int ii=1;ii<Last_x;ii++)
  {

    sm_psi[nii][tjj][tkk] = psi(ii+1,tjj,tkk);
    __syncthreads();

    // contribution to (ii-1)
    // use link loaded previous
    // y face current
    // tjj=0 calc face at 0-1 and tjj=30 calc face at 30-31
  
    if ((tkk>0) && (tkk<nLast_z) && (tjj<nLast_y))
    {
      Real xd=(-V0(tjj,tkk) - V0(1 + tjj,tkk) + V2(tjj,tkk) + V2(1 + tjj,tkk))/4.;
      Real yd=-V1(tjj,tkk) + V1(1 + tjj,tkk);
      Real zd=(-V1(tjj,-1 + tkk) + V1(tjj,1 + tkk) - V1(1 + tjj,-1 + tkk) + V1(1 + tjj,1 + tkk))/4.;


      ycharge = sigmaY(ii,tjj+1,tkk,0) * xd + sigmaY(ii,tjj+1,tkk,1) * yd + sigmaY(ii,tjj+1,tkk,2) * zd ; 


      dV += ycharge;
      sm_psi[3][tjj][tkk]=ycharge;
    }
    __syncthreads();

    if ((tkk>0) && (tkk<nLast_z) && (tjj>0) && (tjj<nLast_y))
      dV -= sm_psi[3][tjj-1][tkk];  //bring from left

    __syncthreads();

    // z face current
    // tkk=0 calc face at 0-1 and tkk=30 calc face at 30-31
    if ((tkk<nLast_z) && (tjj>0) && (tjj<nLast_y))
    {

      Real xd=(-V0(tjj,tkk) - V0(tjj,1 + tkk) + V2(tjj,tkk) + V2(tjj,1 + tkk))/4.;
      Real yd=(-V1(-1 + tjj,tkk) - V1(-1 + tjj,1 + tkk) + V1(1 + tjj,tkk) + V1(1 + tjj,1 + tkk))/4.;
      Real zd=-V1(tjj,tkk) + V1(tjj,1 + tkk);

      zcharge = sigmaZ(ii,tjj,tkk+1,0) * xd + sigmaZ(ii,tjj,tkk+1,1) * yd + sigmaZ(ii,tjj,tkk+1,2) * zd ; 
      dV += zcharge;
      sm_psi[3][tjj][tkk]=zcharge;
    }

    __syncthreads();

    if ((tkk>0) && (tkk<nLast_z) && (tjj>0) && (tjj<nLast_y))
      dV -= sm_psi[3][tjj][tkk-1];

    //__syncthreads();

    // x face current
    if ((tkk>0) && (tkk<nLast_z) && (tjj>0) && (tjj<nLast_y))
    {
      Real xd=-V1(tjj,tkk) + V2(tjj,tkk);
      Real yd=(-V1(-1 + tjj,tkk) + V1(1 + tjj,tkk) - V2(-1 + tjj,tkk) + V2(1 + tjj,tkk))/4.;
      Real zd=(-V1(tjj,-1 + tkk) + V1(tjj,1 + tkk) - V2(tjj,-1 + tkk) + V2(tjj,1 + tkk))/4.;
 
      xcharge = sigmaX(ii+1,tjj,tkk,0) * xd + sigmaX(ii+1,tjj,tkk,1) * yd + sigmaX(ii+1,tjj,tkk,2) * zd ; 
      dV += xcharge;
      //store dV
      npsi(ii,tjj,tkk) = dV;

      dV = -xcharge; //pass to the next cell in x-dir
    }
    tii=pii; pii=cii; cii=nii; nii=tii;
  }
//  #undef V0(y,z)
//  #undef V1(y,z)
//  #undef V2(y,z)
//  #undef sigmaX(x,y,z,dir) 
//  #undef sigmaY(x,y,z,dir) 
//  #undef sigmaZ(x,y,z,dir) 
//  #undef psi(x,y,z) 
//  #undef npsi(x,y,z) 
}

__global__ void map_dVm(wo_array_ptr<Real> dVmOut, ro_array_ptr<Real> dVmRaw, ro_array_ptr<int> remap,int nLocal)
{
  int idx0 = threadIdx.x + blockDim.x*blockIdx.x;
  int stride = blockDim.x * gridDim.x;
  for(int idx = idx0 ; idx<nLocal ; idx+=stride)
  {
      dVmOut[idx] = dVmRaw[remap[idx]];
  }

}

//__global__ void map_V(double * VT, double* V, const int *remap,int nCells)
//{
//  int idx0 = threadIdx.x + blockDim.x*blockIdx.x;
//  int stride = blockDim.x * gridDim.x;
//  for(int idx = idx0 ; idx<nCells ; idx+=stride)
//      VT[remap[idx]] = V[idx];
//}

void call_cuda_kernels(ro_mgarray_ptr<Real> Vm, rw_mgarray_ptr<Real> dVm, ro_mgarray_ptr<Real> sigma, int nx, int ny, int nz, wo_mgarray_ptr<Real> dVmOut, ro_mgarray_ptr<int> lookup,int nLocal)
{
   //determine block dim
   //1. blockdim.z and blockdim.y are determined in a simple way.
   int bdimz = (int)((nz-2)/30) + ((nz-2)%30==0?0:1);
   int bdimy = (int)((ny-2)/30) + ((ny-2)%30==0?0:1);
   int bdimx = (int)((nx-2)/XTILE) + ((nx-2)%XTILE==0?0:1);
   
#ifdef GPU_SM_70
   cudaFuncSetAttribute(diff_6face_v1, cudaFuncAttributePreferredSharedMemoryCarveout, 50);
#endif
   //map_V<<<112,512>>>(VmBlockRaw,VmRaw,lookup,nCells);
   diff_6face_v1<<<dim3(bdimx,bdimy,bdimz),dim3(32,32,1)>>>(&Vm.useOn(GPU)[0],&dVm.useOn(GPU)[0],&sigma.useOn(GPU)[0],&sigma.useOn(GPU)[0]+3*nx*ny*nz,&sigma.useOn(GPU)[0]+6*nx*ny*nz,nx,ny,nz);

//   cudaMemcpy(tmp,dVmRaw,nx*ny*nz*sizeof(double), cudaMemcpyDeviceToHost);
//   for(int ii=0;ii<nx;ii++)
//   for(int jj=0;jj<ny;jj++)
//   for(int kk=0;kk<nz;kk++)
//   {
//     printf("mcp (%d,%d,%d)=%e\n",ii,jj,kk,tmp[kk+nz*(jj+ny*ii)]);
//   }

   map_dVm<<<112,512>>>(dVmOut.useOn(GPU),dVm.useOn(GPU),lookup.useOn(GPU),nLocal);
}

__global__ void copy_to_block_kernel(wo_array_ptr<double> block, ro_array_ptr<int> lookup, ro_array_ptr<double> source, const int begin, const int end)
{
   int ii = threadIdx.x + blockIdx.x*blockDim.x;
   if (ii >= end-begin) { return; }
   block[lookup[ii+begin]] = source[ii];
}

void copy_to_block(wo_mgarray_ptr<double> block, ro_mgarray_ptr<int> lookup, ro_mgarray_ptr<double> source, const int begin, const int end)
{
   int blockSize = 1024;

   if(source.size()>0){
   	copy_to_block_kernel<<<(end-begin+blockSize-1)/blockSize,blockSize>>>
      		(block.useOn(GPU),
       		lookup.useOn(GPU),
       		source.useOn(GPU),
       		begin,
       		end);

    }
}
