//#include "CUDADiffusion.hh"
//#include "DiffusionUtils.hh"
//#include "SymmetricTensor.hh"
//#include <vector>
//#include <map>

//#include "options.h"
//#include "cudautil.h"

typedef double Real;


__global__ void diff_6face_v1(const Real* d_psi, Real* d_npsi, const Real* d_sigmaX, const Real* d_sigmaY, const Real* d_sigmaZ,int Lii, int Ljj, int Lkk)
{

  //rough sizing may be needed
  //e.g. 512x512x512 per GPU
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

  #define psi(x,y,z) d_psi[ z + Lkk * ( y + Ljj * x ) ]
  #define npsi(x,y,z) d_npsi[ z + Lkk * ( y + Ljj * x ) ]

  int tjj = threadIdx.y;
  int tkk = threadIdx.x;

  //shift for each tile
//  d_psi    += 30 * blockIdx.x + Lkk * ( 30 * blockIdx.y );
//  d_npsi   += 30 * blockIdx.x + Lkk * ( 30 * blockIdx.y );
  d_psi    = &(psi(30*blockIdx.x, 30*blockIdx.y, 30*blockIdx.z));
  d_npsi   = &(npsi(30*blockIdx.x, 30*blockIdx.y, 30*blockIdx.z));

  d_sigmaX  = &(sigmaX(30*blockIdx.x, 30*blockIdx.y, 30*blockIdx.z, 0));
  d_sigmaY  = &(sigmaY(30*blockIdx.x, 30*blockIdx.y, 30*blockIdx.z, 0));
  d_sigmaZ  = &(sigmaZ(30*blockIdx.x, 30*blockIdx.y, 30*blockIdx.z, 0));

  int Last_x=31;
  if (blockIdx.x == gridDim.x-1) Last_x = Lii-2 - 30 * blockIdx.x + 1;
//  d_sigmaX += 30 * blockIdx.x + (Lkk-2) * ( 31 * blockIdx.y );
//  d_sigmaY += 30 * blockIdx.x + (Lkk-2) * ( 31 * blockIdx.y );
//  d_sigmaZ += 31 * blockIdx.x + (Lkk-1) * ( 31 * blockIdx.y );

//  printf("tjj tkk bx by = %d %d %d %d\n",tjj,tkk,blockIdx.x,blockIdx.y);


  int pii,cii,nii,tii;
  pii=0; cii=1; nii=2;

  sm_psi[cii][tkk][tjj] = psi(0,tjj,tkk);
  sm_psi[nii][tkk][tjj] = psi(1,tjj,tkk);
  Real xcharge,ycharge,zcharge,dV;

  __syncthreads();
  //initial
  if ((tkk>0) && (tkk<31) && (tjj>0) && (tjj<31))
  {
    Real xd=-V1(tjj,tkk) + V2(tjj,tkk);
    Real yd=(-V1(-1 + tjj,tkk) + V1(1 + tjj,tkk) - V2(-1 + tjj,tkk) + V2(1 + tjj,tkk))/4.;
    Real zd=(-V1(tjj,-1 + tkk) + V1(tjj,1 + tkk) - V2(tjj,-1 + tkk) + V2(tjj,1 + tkk))/4.;

    dV = sigmaX(0,tjj,tkk,0) * xd + sigmaX(0,tjj,tkk,1) * yd + sigmaX(0,tjj,tkk,2) * zd ; 
  }

  tii=pii; pii=cii; cii=nii; nii=tii;

  for(int ii=1;ii<Last_x;ii++)
  {
    sm_psi[nii][tkk][tjj] = psi(ii+1,tjj,tkk);
    __syncthreads();

    // contribution to (ii-1)
    // use link loaded previous
    // y face current
    // tjj=0 calc face at 0-1 and tjj=30 calc face at 30-31
  
    if ((tkk>0) && (tkk<31) && (tjj<31))
    {
      Real xd=(-V0(tjj,tkk) - V0(1 + tjj,tkk) + V2(tjj,tkk) + V2(1 + tjj,tkk))/4.;
      Real yd=-V1(tjj,tkk) + V1(1 + tjj,tkk);
      Real zd=(-V1(tjj,-1 + tkk) + V1(tjj,1 + tkk) - V1(1 + tjj,-1 + tkk) + V1(1 + tjj,1 + tkk))/4.;

      ycharge = sigmaY(ii,tjj,tkk,0) * xd + sigmaY(ii,tjj,tkk,1) * yd + sigmaY(ii,tjj,tkk,2) * zd ; 
      dV += ycharge;
      sm_psi[3][tjj][tkk]=ycharge;
    }
    __syncthreads();

    if ((tkk>0) && (tkk<31) && (tjj>0) && (tjj<31))
      dV -= sm_psi[3][tjj-1][tkk];

    __syncthreads();

    // z face current
    // tkk=0 calc face at 0-1 and tkk=30 calc face at 30-31
    if ((tkk<31) && (tjj>0) && (tjj<31))
    {

      Real xd=(-V0(tjj,tkk) - V0(tjj,1 + tkk) + V2(tjj,tkk) + V2(tjj,1 + tkk))/4.;
      Real yd=(-V1(-1 + tjj,tkk) - V1(-1 + tjj,1 + tkk) + V1(1 + tjj,tkk) + V1(1 + tjj,1 + tkk))/4.;
      Real zd=-V1(tjj,tkk) + V1(tjj,1 + tkk);

      zcharge = sigmaZ(ii,tjj,tkk,0) * xd + sigmaZ(ii,tjj,tkk,1) * yd + sigmaZ(ii,tjj,tkk,2) * zd ; 
      dV += zcharge;
      sm_psi[3][tjj][tkk]=zcharge;
    }

    __syncthreads();

    if ((tkk>0) && (tkk<31) && (tjj>0) && (tjj<30))
      dV -= sm_psi[3][tjj][tkk-1];

    //__syncthreads();

    // x face current
    if ((tkk>0) && (tkk<31) && (tjj>0) && (tjj<31))
    {
      Real xd=-V1(tjj,tkk) + V2(tjj,tkk);
      Real yd=(-V1(-1 + tjj,tkk) + V1(1 + tjj,tkk) - V2(-1 + tjj,tkk) + V2(1 + tjj,tkk))/4.;
      Real zd=(-V1(tjj,-1 + tkk) + V1(tjj,1 + tkk) - V2(tjj,-1 + tkk) + V2(tjj,1 + tkk))/4.;
 
      xcharge = sigmaX(ii,tjj,tkk,0) * xd + sigmaX(ii,tjj,tkk,1) * yd + sigmaX(ii,tjj,tkk,2) * zd ; 
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

__global__ void map_dVm(double * dVmT, double* dVm, const int *remap,int nCells)
{
  int idx0 = threadIdx.x + blockDim.x*blockIdx.x;
  int stride = blockDim.x * gridDim.x;
  for(int idx = idx0 ; idx<nCells ; idx+=stride)
      dVmT[idx] = dVm[remap[idx]];
}


extern "C"
{
void call_cuda_kernels(const Real *VmRaw, Real *dVmRaw, const Real *sigmaRaw, int nx, int ny, int nz, Real *dVmOut, const int *lookup,int nCells)
{
   diff_6face_v1<<<dim3(10,10,10),dim3(32,32,1)>>>(VmRaw,dVmRaw,sigmaRaw,sigmaRaw+3*nx*ny*nz,sigmaRaw+6*nx*ny*nz,nx,ny,nz);
   map_dVm<<<112,512>>>(dVmRaw,dVmOut,lookup,nCells);
}
}
