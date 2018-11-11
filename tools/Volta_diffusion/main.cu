#include <stdio.h>

#define XTILE 20
typedef double Real;

extern "C"
{
void call_cuda_kernels(const Real *VmRaw, Real *dVmRaw, const Real *sigmaRaw, int nx, int ny, int nz, Real *dVmOut, const int *lookup,int nCells);
}

int main()
{
  int nx=100;
  int ny=100;
  int nz=100;

  Real *d_Vm,*d_dVm,*sigmaRaw,*dVmOut;
  Real *h_sigma, *h_Vm, *h_VmRaw;

  cudaMalloc(&d_Vm,sizeof(double)*nx*ny*nz);
  h_Vm=(Real*)malloc(sizeof(double)*nx*ny*nz);

  #define h_Vm(x,y,z) h_Vm[ z + nz * ( y + ny * ( x  ) ) ]

  for(int ii=0;ii<nx;ii++)
  for(int jj=0;jj<ny;jj++)
  for(int kk=0;kk<nz;kk++)
  {
    h_Vm(ii,jj,kk) = sin(ii+jj+kk);
  }

  cudaMemcpy( d_Vm, h_Vm, sizeof(double) * nx*ny*nz , cudaMemcpyHostToDevice );

  cudaMalloc(&sigmaRaw,sizeof(double)*nx*ny*nz*9);
  cudaMemset(sigmaRaw,1,sizeof(double)*nx*ny*nz*9);

  cudaMalloc(&d_dVm,sizeof(double)*nx*ny*nz);
  cudaMemset(d_dVm,1,sizeof(double)*nx*ny*nz);

  cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);

  printf("running with nx=%d ny=%d nz=%d\n",nx,ny,nz);
  call_cuda_kernels(d_Vm,d_dVm,sigmaRaw,nx,ny,nz,0,0,nx*ny*nz);


//  h_sigma=(Real*)malloc(sizeof(double)*nx*ny*nz*9);
//  h_dVm=(Real*)malloc(sizeof(double)*nx*ny*nz);
//  cudaMemset(VmRaw,0,sizeof(double)*nx*ny*nz);
//  cudaMemset(dVmRaw,0,sizeof(double)*nx*ny*nz);
//  cudaMemset(sigmaRaw,0,sizeof(double)*nx*ny*nz*9);
//
//  //test 1
//  //set 1 at ii,jj,kk
//  //set sigmaXYZ 1
//  Real value=1.;
//  int ii=jj=kk=3;
//  cudaMemcpy( &(VmRaw[ kk + nz*( jj + ny* ii ) ]) , &value, sizeof(double), , cudaMemcpyHostToDevice );
//
//  #define sigmaX(x,y,z,dir) h_sigma[ z + Lkk * ( y + Ljj * ( x + Lii * dir ) ) ]
//  #define sigmaY(x,y,z,dir) h_sigma[ z + Lkk * ( y + Ljj * ( x + Lii * dir ) ) + 3*nx*ny*nz ]
//  #define sigmaZ(x,y,z,dir) h_sigma[ z + Lkk * ( y + Ljj * ( x + Lii * dir ) ) + 6*nx*ny*nz ]
//
//  for(int ii=0;ii<nx;ii++)
//  for(int jj=0;jj<ny;jj++)
//  for(int kk=0;kk<nz;kk++)
//  {
//      sigmaX(ii,jj,kk,0) = 1;
//      sigmaY(ii,jj,kk,1) = 1;
//      sigmaZ(ii,jj,kk,2) = 1;
//  }
//      
//  cudaMemcpy( sigmaRaw, h_sigma, sizeof(double) * 9 * nx*ny*nz , cudaMemcpyHostToDevice );
//
//  call_cuda_kernels(VmRaw,dVmRaw,sigmaRaw,nx,ny,nz,0,0,nx*ny*nz);
//
//  cudaMemcpy( h_dVm, dVmRaw, sizeof(double)  * nx*ny*nz , cudaMemcpyDeviceToHost );
//
//  for(int ii=2;ii<5;ii++)
//  for(int jj=2;jj<5;jj++)
//  for(int kk=2;kk<5;kk++)
//  {
//     printf("(%d,%d,%d) = %f\n",ii,jj,kk,h_dVm[ kk + nz*(jj + ny*ii)]);
//  }
  
  cudaDeviceSynchronize();
  cudaFree(d_Vm);
  cudaFree(d_dVm);
  cudaFree(sigmaRaw);
//  free(h_dVm);
  free(h_Vm);
}
