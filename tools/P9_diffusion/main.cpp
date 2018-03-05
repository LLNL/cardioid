#include <sched.h>
// -*-c++-*-

#include <iostream>
#include <vector>
#include "options.h"
#include "types.h"
#include <sys/time.h>
#include <omp.h>

extern "C" {
void HPM_Init(void);
void HPM_Start(const char *);
void HPM_Stop(const char *);
void HPM_Print(void);
}

using namespace std;
void nosimd_diff_cpu(Real* d_psi, Real* d_npsi, Real* d_sigmaX, Real* d_sigmaY, Real* d_sigmaZ,int Lii, int Ljj, int Lkk);
void simd_baseline(Real* psi, Real* npsi, Real* sigmaX, Real* sigmaY, Real* sigmaZ,int Lii, int Ljj, int Lkk);
void simd_hoist(Real* psi, Real* npsi, Real* sigmaX, Real* sigmaY, Real* sigmaZ,int Lii, int Ljj, int Lkk);
void simd_hoist_to_thr(Real* psi, Real* npsi, Real* sigmaX, Real* sigmaY, Real* sigmaZ,int Lii, int Ljj, int Lkk);
void simd_diff_cpu_thr_30x30(Real* psi, Real* npsi, Real* sigmaX, Real* sigmaY, Real* sigmaZ,int Lii, int Ljj, int Lkk);
void simd_mthrd_30disj(int nx,int ny,int nz,int numtimes);

void simd_diff_cpu_thr(Real* psi, Real* npsi, Real* sigmaX, Real* sigmaY, Real* sigmaZ,int Lii, int Ljj, int Lkk,int);
void simd_diff_cpu_thr_lg(Real* psi, Real* npsi, Real* sigmaX, Real* sigmaY, Real* sigmaZ,int Lii, int Ljj, int Lkk,int);
void simd_diff_cpu_thr_30x30_oh(Real* psi, Real* npsi, Real* sigmaX, Real* sigmaY, Real* sigmaZ,int Lii, int Ljj, int Lkk,int);
void simd_diff_cpu_thr_disj(int nx,int ny,int nz,int numtimes);


int main(int argc, char* argv[]) {

   struct gengetopt_args_info params;
   cmdline_parser(argc, argv, &params);

   int nx = params.nx_arg;
   int ny = params.ny_arg;
   int nz = params.nz_arg;

   string kernel = params.kernel_arg;
   int numpoints = nx*ny*nz;
   Real* c_Vm;
   Real* c_dVm;
   
   int numtimes = params.numtimes_arg;
   
   //all kernels should only run on the interior.  Ignore boundary for now.
   if (kernel == "simd_cpu")
   {
     Real* Vm = (Real*)malloc(nx*ny*nz*sizeof(Real));
     Real* dVm = (Real*)malloc(nx*ny*nz*sizeof(Real));

     Real* sx = (Real*)malloc(nx*ny*nz*3*sizeof(Real));
     Real* sy = (Real*)malloc(nx*ny*nz*3*sizeof(Real));
     Real* sz = (Real*)malloc(nx*ny*nz*3*sizeof(Real));
     
     cout << "initializing Vm,dVM,sx,sy,sz" << endl;
     for(int ii=0;ii<nx*ny*nz;ii++) Vm[ii] = 0.21*(ii%10);
     for(int ii=0;ii<nx*ny*nz;ii++) dVm[ii] = 0;
     for(int ii=0;ii<3*nx*ny*nz;ii++) sx[ii] = 0.3*(ii%4);
     for(int ii=0;ii<3*nx*ny*nz;ii++) sy[ii] = 0.4*(ii%4);
     for(int ii=0;ii<3*nx*ny*nz;ii++) sz[ii] = 0.2*(ii%4);

     cout << "simdized 6face on cpu: "<< (nx-2)*(ny-2)*(nz-2) << " cells. " << endl;

     
     struct timeval t0;
     struct timeval t1;

   //HPM_Init();
   //HPM_Start("loop");
     gettimeofday(&t0,NULL);
     for (int ii=0; ii<numtimes; ii++) {
       simd_baseline(Vm,dVm,sx,sy,sz,nx,ny,nz);
     }
     gettimeofday(&t1,NULL);
   //HPM_Stop("loop");
   //HPM_Print();

     cout << "done in " << (t1.tv_sec-t0.tv_sec)*1000000 + t1.tv_usec-t0.tv_usec << " us sec" << endl;

//     cout << "dVm=";
//     for(int ii=0;ii<10;ii++)
//     for(int jj=0;jj<10;jj++)
//     for(int kk=0;kk<10;kk++)
//     {
//       cout << dVm[ii + nx*(jj + ny * kk)] << " "; 
//     }
//     cout << endl;

     free(Vm);
     free(dVm);
     free(sx);
     free(sy);
     free(sz);

   }
   else if (kernel == "simd_cpu_hoist")
   {
     Real* Vm = (Real*)malloc(nx*ny*nz*sizeof(Real));
     Real* dVm = (Real*)malloc(nx*ny*nz*sizeof(Real));

     Real* sx = (Real*)malloc(nx*ny*nz*3*sizeof(Real));
     Real* sy = (Real*)malloc(nx*ny*nz*3*sizeof(Real));
     Real* sz = (Real*)malloc(nx*ny*nz*3*sizeof(Real));
     
     cout << "initializing Vm,dVM,sx,sy,sz" << endl;
     for(int ii=0;ii<nx*ny*nz;ii++) Vm[ii] = 0.21*(ii%10);
     for(int ii=0;ii<nx*ny*nz;ii++) dVm[ii] = 0;
     for(int ii=0;ii<3*nx*ny*nz;ii++) sx[ii] = 0.3*(ii%4);
     for(int ii=0;ii<3*nx*ny*nz;ii++) sy[ii] = 0.4*(ii%4);
     for(int ii=0;ii<3*nx*ny*nz;ii++) sz[ii] = 0.2*(ii%4);

     cout << "simdized hoist on cpu: "<< (nx-2)*(ny-2)*(nz-2) << " cells. " << endl;

     
     struct timeval t0;
     struct timeval t1;

//   HPM_Init();
//   HPM_Start("loop");
     gettimeofday(&t0,NULL);
     for (int ii=0; ii<numtimes; ii++) {
       simd_hoist(Vm,dVm,sx,sy,sz,nx,ny,nz);
     }
     gettimeofday(&t1,NULL);
//   HPM_Stop("loop");
//   HPM_Print();

     cout << "done in " << (t1.tv_sec-t0.tv_sec)*1000000 + t1.tv_usec-t0.tv_usec << " us sec" << endl;

//     cout << "dVm=";
//     for(int ii=0;ii<10;ii++)
//     for(int jj=0;jj<10;jj++)
//     for(int kk=0;kk<10;kk++)
//     {
//       cout << dVm[ii + nx*(jj + ny * kk)] << " "; 
//     }
//     cout << endl;

     free(Vm);
     free(dVm);
     free(sx);
     free(sy);
     free(sz);

   }
   else if (kernel == "simd_hoist_to_mthrd")
   {
     Real* Vm = (Real*)malloc(nx*ny*nz*sizeof(Real));
     Real* dVm = (Real*)malloc(nx*ny*nz*sizeof(Real));

     Real* sx = (Real*)malloc(nx*ny*nz*3*sizeof(Real));
     Real* sy = (Real*)malloc(nx*ny*nz*3*sizeof(Real));
     Real* sz = (Real*)malloc(nx*ny*nz*3*sizeof(Real));
     
     cout << "initializing Vm,dVM,sx,sy,sz" << endl;
     for(int ii=0;ii<nx*ny*nz;ii++) Vm[ii] = 0.21*(ii%10);
     for(int ii=0;ii<nx*ny*nz;ii++) dVm[ii] = 0;
     for(int ii=0;ii<3*nx*ny*nz;ii++) sx[ii] = 0.3*(ii%4);
     for(int ii=0;ii<3*nx*ny*nz;ii++) sy[ii] = 0.4*(ii%4);
     for(int ii=0;ii<3*nx*ny*nz;ii++) sz[ii] = 0.2*(ii%4);

     cout << "simdized hoist on cpu: "<< (nx-2)*(ny-2)*(nz-2) << " cells. " << endl;

     
     struct timeval t0;
     struct timeval t1;

   //HPM_Init();
   //HPM_Start("loop");
     simd_hoist_to_thr(Vm,dVm,sx,sy,sz,nx,ny,nz);
     gettimeofday(&t0,NULL);
     for (int ii=0; ii<numtimes; ii++) {
       simd_hoist_to_thr(Vm,dVm,sx,sy,sz,nx,ny,nz);
     }
     gettimeofday(&t1,NULL);
   //HPM_Stop("loop");
   //HPM_Print();

     cout << "done in " << (t1.tv_sec-t0.tv_sec)*1000000 + t1.tv_usec-t0.tv_usec << " us sec" << endl;

//     cout << "dVm=";
//     for(int ii=0;ii<10;ii++)
//     for(int jj=0;jj<10;jj++)
//     for(int kk=0;kk<10;kk++)
//     {
//       cout << dVm[ii + nx*(jj + ny * kk)] << " "; 
//     }
//     cout << endl;

     free(Vm);
     free(dVm);
     free(sx);
     free(sy);
     free(sz);

   }
   else if (kernel == "simd_cpu_thr")
   {
     Real* Vm = (Real*)malloc(nx*ny*nz*sizeof(Real));
     Real* dVm = (Real*)malloc(nx*ny*nz*sizeof(Real));

     Real* sx = (Real*)malloc(nx*ny*nz*3*sizeof(Real));
     Real* sy = (Real*)malloc(nx*ny*nz*3*sizeof(Real));
     Real* sz = (Real*)malloc(nx*ny*nz*3*sizeof(Real));
     
     cout << "initializing Vm,dVM,sx,sy,sz" << endl;
     for(int ii=0;ii<nx*ny*nz;ii++) Vm[ii] = 0.21*(ii%10);
     for(int ii=0;ii<nx*ny*nz;ii++) dVm[ii] = 0;
     for(int ii=0;ii<3*nx*ny*nz;ii++) sx[ii] = 0.3*(ii%4);
     for(int ii=0;ii<3*nx*ny*nz;ii++) sy[ii] = 0.4*(ii%4);
     for(int ii=0;ii<3*nx*ny*nz;ii++) sz[ii] = 0.2*(ii%4);

     cout << "simdized threaded on cpu: "<< (nx-2)*(ny-2)*(nz-2) << " cells. " << endl;
     struct timeval t0;
     struct timeval t1;

   //HPM_Init();
   //HPM_Start("loop");
     gettimeofday(&t0,NULL);
     simd_diff_cpu_thr(Vm,dVm,sx,sy,sz,nx,ny,nz,numtimes);
     gettimeofday(&t1,NULL);
   //HPM_Stop("loop");
   //HPM_Print();
     

     cout << "done in " << (t1.tv_sec-t0.tv_sec)*1000000 + t1.tv_usec-t0.tv_usec << " us sec" << endl;

     free(Vm);
     free(dVm);
     free(sx);
     free(sy);
     free(sz);

   }
   else if (kernel == "simd_cpu_thr_30")
   {
     Real* Vm = (Real*)malloc(nx*ny*nz*sizeof(Real));
     Real* dVm = (Real*)malloc(nx*ny*nz*sizeof(Real));

     Real* sx = (Real*)malloc(nx*ny*nz*3*sizeof(Real));
     Real* sy = (Real*)malloc(nx*ny*nz*3*sizeof(Real));
     Real* sz = (Real*)malloc(nx*ny*nz*3*sizeof(Real));
     
     cout << "initializing Vm,dVM,sx,sy,sz" << endl;
     for(int ii=0;ii<nx*ny*nz;ii++) Vm[ii] = 0.21*(ii%10);
     for(int ii=0;ii<nx*ny*nz;ii++) dVm[ii] = 0;
     for(int ii=0;ii<3*nx*ny*nz;ii++) sx[ii] = 0.3*(ii%4);
     for(int ii=0;ii<3*nx*ny*nz;ii++) sy[ii] = 0.4*(ii%4);
     for(int ii=0;ii<3*nx*ny*nz;ii++) sz[ii] = 0.2*(ii%4);

     cout << "simdized 30x30 on cpu: "<< (nx-2)*(ny-2)*(nz-2) << " cells. " << endl;
     struct timeval t0;
     struct timeval t1;

   //HPM_Init();
   //HPM_Start("loop");
     simd_diff_cpu_thr_30x30(Vm,dVm,sx,sy,sz,nx,ny,nz);
     gettimeofday(&t0,NULL);
     for (int ii=0; ii<numtimes; ii++) {
       simd_diff_cpu_thr_30x30(Vm,dVm,sx,sy,sz,nx,ny,nz);
     }
     gettimeofday(&t1,NULL);
   //HPM_Stop("loop");
   //HPM_Print();
     

     cout << "done in " << (t1.tv_sec-t0.tv_sec)*1000000 + t1.tv_usec-t0.tv_usec << " us sec" << endl;

//     cout << "dVm=";
//     for(int ii=0;ii<10;ii++)
//     for(int jj=0;jj<10;jj++)
//     for(int kk=0;kk<10;kk++)
//     {
//       cout << dVm[ii + nx*(jj + ny * kk)] << " "; 
//     }

     cout << endl;
     free(Vm);
     free(dVm);
     free(sx);
     free(sy);
     free(sz);

   }
   else if (kernel == "simd_cpu_thr_30_oh")
   {
     Real* Vm = (Real*)malloc(nx*ny*nz*sizeof(Real));
     Real* dVm = (Real*)malloc(nx*ny*nz*sizeof(Real));

     Real* sx = (Real*)malloc(nx*ny*nz*3*sizeof(Real));
     Real* sy = (Real*)malloc(nx*ny*nz*3*sizeof(Real));
     Real* sz = (Real*)malloc(nx*ny*nz*3*sizeof(Real));
     
     cout << "initializing Vm,dVM,sx,sy,sz" << endl;
     for(int ii=0;ii<nx*ny*nz;ii++) Vm[ii] = 0.21*(ii%10);
     for(int ii=0;ii<nx*ny*nz;ii++) dVm[ii] = 0;
     for(int ii=0;ii<3*nx*ny*nz;ii++) sx[ii] = 0.3*(ii%4);
     for(int ii=0;ii<3*nx*ny*nz;ii++) sy[ii] = 0.4*(ii%4);
     for(int ii=0;ii<3*nx*ny*nz;ii++) sz[ii] = 0.2*(ii%4);

     cout << "simdized 30x30 on cpu: "<< (nx-2)*(ny-2)*(nz-2) << " cells. " << endl;
     struct timeval t0;
     struct timeval t1;

   //HPM_Init();
   //HPM_Start("loop");
     gettimeofday(&t0,NULL);
     simd_diff_cpu_thr_30x30_oh(Vm,dVm,sx,sy,sz,nx,ny,nz,numtimes);
     gettimeofday(&t1,NULL);
   //HPM_Stop("loop");
   //HPM_Print();
     

     cout << "done in " << (t1.tv_sec-t0.tv_sec)*1000000 + t1.tv_usec-t0.tv_usec << " us sec" << endl;

//     cout << "dVm=";
//     for(int ii=0;ii<10;ii++)
//     for(int jj=0;jj<10;jj++)
//     for(int kk=0;kk<10;kk++)
//     {
//       cout << dVm[ii + nx*(jj + ny * kk)] << " "; 
//     }

     cout << endl;
     free(Vm);
     free(dVm);
     free(sx);
     free(sy);
     free(sz);

   }
   else if (kernel == "simd_cpu_thr_lg")
   {
     Real* Vm = (Real*)malloc(nx*ny*nz*sizeof(Real));
     Real* dVm = (Real*)malloc(nx*ny*nz*sizeof(Real));

     Real* sx = (Real*)malloc(nx*ny*nz*3*sizeof(Real));
     Real* sy = (Real*)malloc(nx*ny*nz*3*sizeof(Real));
     Real* sz = (Real*)malloc(nx*ny*nz*3*sizeof(Real));
     
     cout << "initializing Vm,dVM,sx,sy,sz" << endl;
     for(int ii=0;ii<nx*ny*nz;ii++) Vm[ii] = 0.21*(ii%10);
     for(int ii=0;ii<nx*ny*nz;ii++) dVm[ii] = 0;
     for(int ii=0;ii<3*nx*ny*nz;ii++) sx[ii] = 0.3*(ii%4);
     for(int ii=0;ii<3*nx*ny*nz;ii++) sy[ii] = 0.4*(ii%4);
     for(int ii=0;ii<3*nx*ny*nz;ii++) sz[ii] = 0.2*(ii%4);

     cout << "simd threadedlarge on cpu: "<< (nx-2)*(ny-2)*(nz-2) << " cells. " << endl;
     struct timeval t0;
     struct timeval t1;

   //HPM_Init();
   //HPM_Start("loop");
     gettimeofday(&t0,NULL);
     simd_diff_cpu_thr_lg(Vm,dVm,sx,sy,sz,nx,ny,nz,numtimes);
     gettimeofday(&t1,NULL);
   //HPM_Stop("loop");
   //HPM_Print();
     

     cout << "done in " << (t1.tv_sec-t0.tv_sec)*1000000 + t1.tv_usec-t0.tv_usec << " us sec" << endl;

     free(Vm);
     free(dVm);
     free(sx);
     free(sy);
     free(sz);

   }
   else if (kernel == "simd_mthrd_30disj")
   {
     cout << "simdized 30dijointthread on cpu: "<< (nx-2)*(ny-2)*(nz-2) << " cells. " << endl;
   //HPM_Init();
   //HPM_Start("loop");
     simd_mthrd_30disj(nx,ny,nz,numtimes);

   //HPM_Stop("loop");
   //HPM_Print();
     
   }
   else if (kernel == "simd_cpu_thr2")
   {



     cout << "simdized dijointthread on cpu: "<< (nx-2)*(ny-2)*(nz-2) << " cells. " << endl;
     simd_diff_cpu_thr_disj(nx,ny,nz,numtimes);

   //HPM_Init();
   //HPM_Start("loop");
   //HPM_Stop("loop");
   //HPM_Print();
     
   }
   else
   {
      cout << "Didn't recognize diffusion kernel \"" << kernel << "\"\n"; 
      return EXIT_FAILURE;
   }

   return 0;
}

