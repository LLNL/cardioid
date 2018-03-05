#include <sched.h>
#include "types.h"
#include "options.h"
#include "altivec.h"
//#include <stdlib.h>
//#include <stdio.h>
//#include "veclib.1.0.4/include/veclib_types.h"
//#include "veclib.1.0.4/include/vec256dp.h"
#define load4(x) (*x)
#define store4(x,y) *x=y
#define make4(x) x

typedef Real vector4;

void nosimd_diff_cpu(Real* d_psi, Real* d_npsi, Real* d_sigmaX, Real* d_sigmaY, Real* d_sigmaZ,int Lii, int Ljj, int Lkk)
{

  //rough sizing may be needed
  //map z dir to SIMD
  //z is the fastest varying direction

  //2d decomposition
  //15x15 in y z direction

  #define sigmaX(x,y,z,dir) d_sigmaX[ z + (Lkk-2) * ( y + (Ljj-2) * ( x + (Lii-1) * dir ) ) ]
  #define sigmaY(x,y,z,dir) d_sigmaY[ z + (Lkk-2) * ( y + (Ljj-1) * ( x + (Lii-2) * dir ) ) ]
  #define sigmaZ(x,y,z,dir) d_sigmaZ[ z + (Lkk-1) * ( y + (Ljj-2) * ( x + (Lii-2) * dir ) ) ]

  #define psi(x,y,z) d_psi[ z + Lkk * ( y + Ljj * x ) ]
  #define npsi(x,y,z) d_npsi[ z + Lkk * ( y + Ljj * x ) ]
  #define V(x,y,z) psi(x,y,z)

  vector4 xm[15][16]; //temporary array
  vector4 ym[16]; //temporary array

  if ((Ljj-2)%15 != 0 ) printf("Ljj-2 must be multiple of 15. This will be fixed, later\n");
  if ((Lkk-2)%15 != 0 ) printf("Lkk-2 must be multiple of 15. This will be fixed, later\n");

  //todo: when Ljj is not multiple of 15
  for(int cjj=1;cjj<Ljj-1;cjj+=15)
  for(int ckk=1;ckk<Lkk-1;ckk+=15)
  {

    __memset(xm,0,sizeof(xm));
//    for(int tjj=0;tjj<15;tjj++)
//    for(int tkk=0;tkk<15;tkk++)
//      xm[tjj][tkk] = make4(0);

    for(int x=1;x<Lii-1;x++)
    {
        //z direction: all loads are un-aligned load
        for(int tjj=0;tjj<15;tjj++)       
        {                                 
          vector4 last = 0;

          for(int tkk=0;tkk<16;tkk++)     
          {                               
            int y = cjj + tjj;            
            int z = ckk + tkk;       

            vector4 v00 = load4( &V(x,y,z-1)   );
            vector4 v0p = load4( &V(x,y+1,z-1) );
            vector4 v0m = load4( &V(x,y-1,z-1) );
            vector4 vp0 = load4( &V(x+1,y,z-1) );
            vector4 vm0 = load4( &V(x-1,y,z-1) );
           
            vector4 w00 = load4( &V(x,y,z)   );
            vector4 w0p = load4( &V(x,y+1,z) );
            vector4 w0m = load4( &V(x,y-1,z) );
            vector4 wp0 = load4( &V(x+1,y,z) );
            vector4 wm0 = load4( &V(x-1,y,z) );
            
            vector4 sX = load4( &sigmaZ(x,y,z-1,0) );
            vector4 sY = load4( &sigmaZ(x,y,z-1,1) );
            vector4 sZ = load4( &sigmaZ(x,y,z-1,2) );
           
            vector4 zzm = sZ * (w00 - v00);
            vector4 zxm = sX * ( vp0 - vm0 + wp0 - wm0 ) * make4(0.5);
            vector4 zym = sY * ( v0p - v0m + w0p - w0m ) * make4(0.5);
          
            vector4 tmp = zzm + zzm + zym;
            store4(&npsi(x,y,z), tmp + last);
            last = tmp;
          } 
        } 

        //x direction: all loads are un-aligned load
        for(int tjj=0;tjj<15;tjj++)       
        {                                 
          for(int tkk=0;tkk<16;tkk++)     
          {                               
            int y = cjj + tjj;            
            int z = ckk + tkk;       

            vector4 v0m = load4( &V(x,y,z-1) );
            vector4 v0p = load4( &V(x,y,z+1) );
            vector4 v00 = load4( &V(x,y,z) );

            vector4 vp0 = load4( &V(x,y+1,z) );
            vector4 vm0 = load4( &V(x,y-1,z) );

            vector4 w0m = load4( &V(x+1,y,z-1) );
            vector4 w0p = load4( &V(x+1,y,z+1) );
            vector4 w00 = load4( &V(x+1,y,z) );

            vector4 wp0 = load4( &V(x+1,y+1,z) );
            vector4 wm0 = load4( &V(x+1,y-1,z) );
           
            vector4 sX = load4( &sigmaX(x,y,z,0) );
            vector4 sY = load4( &sigmaX(x,y,z,1) );
            vector4 sZ = load4( &sigmaX(x,y,z,2) );
           
            vector4 xxm = sX * (w00 - v00);
            vector4 xym = sY * (vp0 - vm0 + wp0 - wm0 ) * make4(0.5);
            vector4 xzm = sZ * (v0p - v0m + w0p - w0m ) * make4(0.5);
          
            vector4 tmp = xxm + xym + xzm;
            store4(&npsi(x,y,z), xm[tjj][tkk] + tmp);

            xm[tjj][tkk] = tmp;
          } 
        } 

        //y direction:         
        __memset(ym,0,sizeof(ym));
//        for(int tkk=0;tkk<16;tkk++) ym[tkk]=0;
     
        for(int tjj=0;tjj<15;tjj++)       
        {                                 
          for(int tkk=0;tkk<16;tkk++)     
          {                               
            int y = cjj + tjj;            
            int z = ckk + tkk;       

            vector4 v0m = load4( &V(x,y,z-1) );
            vector4 v0p = load4( &V(x,y,z+1) );
            vector4 v00 = load4( &V(x,y,z) );

            vector4 vp0 = load4( &V(x+1,y,z) );
            vector4 vm0 = load4( &V(x-1,y,z) );

            vector4 w0m = load4( &V(x,y+1,z-1) );
            vector4 w0p = load4( &V(x,y+1,z+1) );
            vector4 w00 = load4( &V(x,y+1,z) );

            vector4 wp0 = load4( &V(x+1,y+1,z) );
            vector4 wm0 = load4( &V(x-1,y+1,z) );
           
            vector4 sX = load4( &sigmaY(x,y,z,0) );
            vector4 sY = load4( &sigmaY(x,y,z,1) );
            vector4 sZ = load4( &sigmaY(x,y,z,2) );
           
            vector4 yxm = sX * (vp0 - vm0 + wp0 - wm0 ) * make4(0.5);
            vector4 yym = sY * (w00 - v00);
            vector4 yzm = sZ * (v0p - v0m + w0p - w0m ) * make4(0.5);
          
            vector4 tmp = yxm + yym + yzm;
            store4(&npsi(x,y,z), ym[tkk] + tmp);
            ym[tkk] = tmp;
          } 
        } 
    } //x loop
  } //cjj,ckk
}
