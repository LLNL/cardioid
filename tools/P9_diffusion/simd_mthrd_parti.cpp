#include <sched.h>
#include "types.h"
#include "options.h"
#include <omp.h>
#include <sys/time.h>
#include <iostream>
#include "altivec.h"
//#include "veclib.1.0.4/include/veclib_types.h"
//#include "veclib.1.0.4/include/vec256dp.h"
#define load4(x) vec_xld2(0,x)
#define store4(x,y) vec_xstd2(y,0,x)
#define make4(x) (vector double)(x)

#define loadV(x,y,z) vec_xld2( 2*(z) + 2*Lkk*(y) + 2*LkkLjj*(x) , baseV )
#define loadS(x,y,z,dir) vec_xld2(2*(z) + 2*(Lkk-1)*(y) + 2*Lkk_1Ljj_1*(x) + 2*Lkk_1Ljj_1Lii_1 * (dir) , baseS)
#define storeV(x,y,z,v) vec_xstd2(v,2*(z) + 2*Lkk*(y) + 2*LkkLjj*(x),baseST)

typedef vector double vector4;

void simd_diff_cpu_thr_disj(int nx, int ny, int nz, int numtimes)
{
  //setting num threads
  if ((ny-2)%15 != 0 ) printf("ny-2 must be multiple of 15. This will be fixed, later\n");
  if ((nz-2)%15 != 0 ) printf("nz-2 must be multiple of 15. This will be fixed, later\n");
  int yblk = (ny-2)/15;  
  int zblk = (nz-2)/15;  
  int tblk = yblk * zblk;

  omp_set_num_threads(tblk);
  printf("using %d threads\n",tblk);

  Real* Vm = (Real*)malloc(tblk*17*17*nx*sizeof(Real));
  Real* dVm = (Real*)malloc(tblk*17*17*nx*sizeof(Real));

  Real* sx = (Real*)malloc(tblk*17*17*nx*3*sizeof(Real));
  Real* sy = (Real*)malloc(tblk*17*17*nx*3*sizeof(Real));
  Real* sz = (Real*)malloc(tblk*17*17*nx*3*sizeof(Real));

  vector4* send_py = (vector4*)malloc(tblk*17*nx*sizeof(Real));
  vector4* send_my = (vector4*)malloc(tblk*17*nx*sizeof(Real));
  vector4* send_pz = (vector4*)malloc(tblk*17*nx*sizeof(Real));
  vector4* send_mz = (vector4*)malloc(tblk*17*nx*sizeof(Real));
    
  struct timeval t0;
  struct timeval t1;

  
  std::cout << "initializing Vm,dVM,sx,sy,sz" << Vm << dVm << sx << sy << sz << std::endl;
  #pragma omp parallel shared(nx) firstprivate(Vm,dVm,sx,sy,sz)
  {
    int tid = omp_get_thread_num();

    Vm +=  tid*17*17*nx;
    dVm += tid*17*17*nx;

    sx += tid*17*17*nx*3;
    sy += tid*17*17*nx*3;
    sz += tid*17*17*nx*3;

    for(int ii=0;ii<nx*17*17;ii++) Vm[ii] = 0.21*(ii%10);
    for(int ii=0;ii<nx*17*17;ii++) dVm[ii] = 0;
    for(int ii=0;ii<3*nx*17*17;ii++) sx[ii] = 0.3*(ii%4);
    for(int ii=0;ii<3*nx*17*17;ii++) sy[ii] = 0.4*(ii%4);
    for(int ii=0;ii<3*nx*17*17;ii++) sz[ii] = 0.2*(ii%4);
  }
  std::cout << "done Vm,dVM,sx,sy,sz" << Vm << dVm << sx << sy << sz << std::endl;

  gettimeofday(&t0,NULL);
  //2d decomposition
  //15x15 in y z direction

  #define sigmaX(t,x,y,z,dir) sx[16*16*nx*(3*(t) + dir) + 2*(z + 16 * (y + 16 * (x) ) )]
  #define sigmaY(t,x,y,z,dir) sy[16*16*nx*(3*(t) + dir) + 2*(z + 16 * (y + 16 * (x) ) )]
  #define sigmaZ(t,x,y,z,dir) sz[16*16*nx*(3*(t) + dir) + 2*(z + 16 * (y + 16 * (x) ) )]

  #define psi(t,x,y,z)   Vm[(t)*17*17*nx + 2*(z + 17 * ( y + 17 * (x) ))]
  #define npsi(t,x,y,z) dVm[(t)*17*17*nx + 2*(z + 17 * ( y + 17 * (x) ))]
  #define V(t,x,y,z) psi(t,x,y,z)

  int Lii = nx/2;
  const int Lkk = nz;
  const int Ljj = ny;

  const int LkkLjj = Lkk * Ljj;
  const int Lkk_1Ljj_1 = (Lkk-1)*(Ljj-1);
  const int Lkk_1Ljj_1Lii_1 = (Lkk-1)*(Ljj-1)*(Lii-1);
 
  int kblocks = (Lkk-2)/15;
  int t_tile = (Ljj-2)*kblocks/15;
     
  //printf("totla number of threads : %d\n", omp_get_num_threads());

  #pragma omp parallel 
  {
    int otn = omp_get_thread_num();
    int ont = omp_get_num_threads();

    const int cjj = 1;
    const int ckk = 1;

    vector4 xm[15][16]; //temporary array
    vector4 ym[16]; //temporary array

    for(int tid=otn;tid<tblk;tid+=ont)
    {
      for (int ii=0; ii<numtimes; ii++)
      {
          __memset(xm,0,sizeof(xm));
      
          for(int x=1;x<Lii-1;x++)
          {
              //z direction: all loads are un-aligned load
              //
      
              for(int tjj=0;tjj<15;tjj++)       
              {                                 
                vector4 last = make4(0);
      
                double *baseV =(double*)(&V(tid,x,cjj+tjj,ckk));
                double *baseS = (double*)(&sigmaZ(tid,x,cjj+tjj,ckk,0));
                double *baseST = (double*)(&npsi(tid,x,cjj+tjj,ckk));
      
                for(int tkk=0;tkk<16;tkk++)     
                {                               
                  int z = tkk;       
      
                  vector4 v00 = loadV(0,0,z-1)   ;
                  vector4 v0p = loadV(0,0+1,z-1) ;
                  vector4 v0m = loadV(0,0-1,z-1) ;
                  vector4 vp0 = loadV(0+1,0,z-1) ;
                  vector4 vm0 = loadV(0-1,0,z-1) ;
                 
                  vector4 w00 = loadV(0,0,z)   ;
                  vector4 w0p = loadV(0,0+1,z) ;
                  vector4 w0m = loadV(0,0-1,z) ;
                  vector4 wp0 = loadV(0+1,0,z) ;
                  vector4 wm0 = loadV(0-1,0,z) ;
                  
                  vector4 sX = loadS( 0,0,z-1,0) ;
                  vector4 sY = loadS( 0,0,z-1,1) ;
                  vector4 sZ = loadS( 0,0,z-1,2) ;
                 
                  vector4 zzm = sZ * (w00 - v00);
                  vector4 zxm = sX * ( vp0 - vm0 + wp0 - wm0 ) * make4(0.5);
                  vector4 zym = sY * ( v0p - v0m + w0p - w0m ) * make4(0.5);
                
                  vector4 tmp = zzm + zzm + zym;
                  storeV(0,0,z, tmp + last);
                  last = tmp;
                } 
              } 
      
              //x direction: all loads are un-aligned load
              for(int tjj=0;tjj<15;tjj++)       
              {                                 
      
                double *baseV =(double*)(&V(tid,x,cjj+tjj,ckk));
                double *baseS = (double*)(&sigmaX(tid,x,cjj+tjj,ckk,0));
                double *baseST = (double*)(&npsi(tid,x,cjj+tjj,ckk));
      
                for(int tkk=0;tkk<16;tkk++)     
                {                               
                  int z = tkk;       
      
                  vector4 v0m = loadV( 0,0,z-1) ;
                  vector4 v0p = loadV( 0,0,z+1) ;
                  vector4 v00 = loadV( 0,0,z) ;
      
                  vector4 vp0 = loadV( 0,0+1,z) ;
                  vector4 vm0 = loadV( 0,0-1,z) ;
      
                  vector4 w0m = loadV( 0+1,0,z-1) ;
                  vector4 w0p = loadV( 0+1,0,z+1) ;
                  vector4 w00 = loadV( 0+1,0,z) ;
      
                  vector4 wp0 = loadV( 0+1,0+1,z) ;
                  vector4 wm0 = loadV( 0+1,0-1,z) ;
                 
                  vector4 sX = loadS( 0,0,z,0) ;
                  vector4 sY = loadS( 0,0,z,1) ;
                  vector4 sZ = loadS( 0,0,z,2) ;
                 
                  vector4 xxm = sX * (w00 - v00);
                  vector4 xym = sY * (vp0 - vm0 + wp0 - wm0 ) * make4(0.5);
                  vector4 xzm = sZ * (v0p - v0m + w0p - w0m ) * make4(0.5);
                
                  vector4 tmp = xxm + xym + xzm;
                  storeV(0,0,z, xm[tjj][tkk] + tmp);
      
                  xm[tjj][tkk] = tmp;
                } 
              } 
      
              //y direction:         
              __memset(ym,0,sizeof(ym));
      //        for(int tkk=0;tkk<16;tkk++) ym[tkk]=make4(0);
           
              for(int tjj=0;tjj<15;tjj++)       
              {                                 
      
                double *baseV =(double*)(&V(tid,x,cjj+tjj,ckk));
                double *baseS = (double*)(&sigmaY(tid,x,cjj+tjj,ckk,0));
                double *baseST = (double*)(&npsi(tid,x,cjj+tjj,ckk));
      
                for(int tkk=0;tkk<16;tkk++)     
                {                               
                  int z = tkk;       
      
                  vector4 v0m = loadV( 0,0,z-1) ;
                  vector4 v0p = loadV( 0,0,z+1) ;
                  vector4 v00 = loadV( 0,0,z) ;
      
                  vector4 vp0 = loadV( 0+1,0,z) ;
                  vector4 vm0 = loadV( 0-1,0,z) ;
      
                  vector4 w0m = loadV( 0,0+1,z-1) ;
                  vector4 w0p = loadV( 0,0+1,z+1) ;
                  vector4 w00 = loadV( 0,0+1,z) ;
      
                  vector4 wp0 = loadV( 0+1,0+1,z) ;
                  vector4 wm0 = loadV( 0-1,0+1,z) ;
                 
                  vector4 sX = loadS (0,0,z,0) ;
                  vector4 sY = loadS (0,0,z,1) ;
                  vector4 sZ = loadS (0,0,z,2) ;
                 
                  vector4 yxm = sX * (vp0 - vm0 + wp0 - wm0 ) * make4(0.5);
                  vector4 yym = sY * (w00 - v00);
                  vector4 yzm = sZ * (v0p - v0m + w0p - w0m ) * make4(0.5);
                
                  vector4 tmp = yxm + yym + yzm;
                  storeV(0,0,z, ym[tkk] + tmp);
                  ym[tkk] = tmp;
                } 
              } 
          } //x loop

          //data prep : pack
          int tt=0;

          for(int x=1;x<Lii-1;x++)
          {
             double *baseV = (double*)(&V(tid,x,0,0));
             for(int tjj=0;tjj<15;tjj++)       
             {                                 
               send_py[tt] = loadV(0,15,tjj);  
               send_my[tt] = loadV(0,1,tjj);  
      
               send_pz[tt] = loadV(0,tjj,15);
               send_mz[tt] = loadV(0,tjj,1);
         
               tt++;
             }
             //4 more to be collected
          }


          //sync and copy
          #pragma omp barrier 
          
          
      } //num times
    } // tid
  } // omp parallel


  gettimeofday(&t1,NULL);
  std::cout << "done in " << (t1.tv_sec-t0.tv_sec)*1000000 + t1.tv_usec-t0.tv_usec << " us sec" << std::endl;

}
