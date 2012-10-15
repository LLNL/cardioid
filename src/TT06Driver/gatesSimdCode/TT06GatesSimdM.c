#include <assert.h>
#include <stdio.h>
#include "portableSIMD.h" 

/*
 * Apparently the 12.1 version of mpixlc_r sometimes improperly implements the vec_lds intrinsic.
 * On some loads, only the first of the 4 doubles of the vector register is actually loaded with the value.
 * As a work-around, we use vec_splats instead.
 * In principle, this is probably slower, since it will presumably first cause the value to be 
 * loaded into a register, then splat'ed into another (or perhaps the same register).
 * Besides involving at least a second instruction, it also creates a delay of several cycles.
 * (Of course, if the compiler were smart, it would be able to convert vec_splats(*mem) into vec_lds(0,mem),
 * but since we get different results, presumably it isn't.
 *
 * At some point, this bug will probably be fixed.  To test find out, just try not defining BAD_VEC_LDS.
 */
#define	BAD_VEC_LDS

#ifdef	BAD_VEC_LDS
// NOTE this will not work for n!=0
#define VEC_LDS(n,addr) vec_splats(*(addr))
#else
#define VEC_LDS vec_lds
#endif

/*
#ifdef BGQ 
#define get
#else
#define get .v
#endif
*/
//#define TTROUND(nc,n) (((nc)-3)*(n)/(n))
//#define TTROUND(nc,n) ((n)*(((nc)-3)/(n)))

#define TTROUND(nc,n) ((nc)-((n)-1))


static  double exp_a[32]__attribute__((aligned(32)));
void initExp()
{
   exp_a[0]=1.0; 
   for (int i=1;i<32;i++) exp_a[i] = exp_a[i-1]/i;
}
void update_mGate_v1(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a)
{
  //update_mGate(dt, nCells, VM, g, mhu_a, tauR_a); return;
typedef vector4double vdt;

int gateIndex=0;
int mhu_l=10;
int mhu_m=5;
int tauR_l=1;
int tauR_m=18;
//int exp_l=16;
int exp_l=12;

/*
  int mhu_k  = 14; 
  int tauR_k = 18; 
*/
   
 double  tauRdt_a[18];

// the number of following v register will be smaller after Jim tune the code
  vdt v_exp_a0  = VEC_LDS(0, &exp_a[0]); 
  vdt v_exp_a1  = VEC_LDS(0, &exp_a[1]);
  vdt v_exp_a2  = VEC_LDS(0, &exp_a[2]); 
  vdt v_exp_a3  = VEC_LDS(0, &exp_a[3]);
  vdt v_exp_a4  = VEC_LDS(0, &exp_a[4]); 
  vdt v_exp_a5  = VEC_LDS(0, &exp_a[5]);
  vdt v_exp_a6  = VEC_LDS(0, &exp_a[6]); 
  vdt v_exp_a7  = VEC_LDS(0, &exp_a[7]);
  vdt v_exp_a8  = VEC_LDS(0, &exp_a[8]); 
  vdt v_exp_a9  = VEC_LDS(0, &exp_a[9]);
  vdt v_exp_a10 = VEC_LDS(0, &exp_a[10]); 
  vdt v_exp_a11 = VEC_LDS(0, &exp_a[11]);
  vdt v_exp_a12 = VEC_LDS(0, &exp_a[12]); 
  //vdt v_exp_a13 = VEC_LDS(0, &exp_a[13]);
  //vdt v_exp_a14 = VEC_LDS(0, &exp_a[14]); 
  //vdt v_exp_a15 = VEC_LDS(0, &exp_a[15]);
  //vdt v_exp_a16 = VEC_LDS(0, &exp_a[16]);


  vdt v_ONE = vec_splats(1.0);

// for (int i=1;i<=exp_l;i++) exp_a[i] = exp_a[i-1]/i;

//  for (int j=17;j>=0;j--)    tauRdt_a[j] = tauR_a[j]*dt;
  tauRdt_a[0]  = tauR_a[0]*dt,
    tauRdt_a[1]  = tauR_a[1]*dt,
   tauRdt_a[2]  = tauR_a[2]*dt,
   tauRdt_a[3]  = tauR_a[3]*dt,
   tauRdt_a[4]  = tauR_a[4]*dt,
   tauRdt_a[5]  = tauR_a[5]*dt,
   tauRdt_a[6]  = tauR_a[6]*dt, 
   tauRdt_a[7]  = tauR_a[7]*dt,
   tauRdt_a[8]  = tauR_a[8]*dt;
   tauRdt_a[9]  = tauR_a[9]*dt,
   tauRdt_a[10] = tauR_a[10]*dt,
   tauRdt_a[11] = tauR_a[11]*dt,
   tauRdt_a[12] = tauR_a[12]*dt,
   tauRdt_a[13] = tauR_a[13]*dt,
   tauRdt_a[14] = tauR_a[14]*dt,
   tauRdt_a[15] = tauR_a[15]*dt, 
   tauRdt_a[16] = tauR_a[16]*dt,
   tauRdt_a[17] = tauR_a[17]*dt;

 int ii = 0;

 for (;ii<TTROUND(nCells,8);ii+=8)
 { 
   vdt v_xa = vec_ld(0, &VM[ii]);  
   vdt v_xb = vec_ld(0, &VM[ii+4]);  

//   for (int j=4;  j>=0; j--)sum1 = mhu_a[j]    + x*sum1;
//   for (int j=14; j>=5; j--)sum2 = mhu_a[j]    + x*sum2;
//   for (int j=17; j>=0; j--)sum3 = tauRdt_a[j] + x*sum3;


//1 
/*
   sum1[ii] = mhu_a[4];
   sum1[ii] = mhu_a[3]     + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[14];
   sum2[ii] = mhu_a[13]    + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[17];
   sum3[ii] = tauRdt_a[16] + VM[ii]*sum3[ii];
*/

   vdt v_mhu_A1    = VEC_LDS(0, &mhu_a[3]);
   vdt v_mhu_A2    = VEC_LDS(0, &mhu_a[4]);
   vdt v_mhu_B1    = VEC_LDS(0, &mhu_a[13]);
   vdt v_mhu_B2    = VEC_LDS(0, &mhu_a[14]);
   vdt v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[16]);
   vdt v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[17]);

   vdt v_sum1a = vec_madd(v_xa, v_mhu_A2,    v_mhu_A1);
   vdt v_sum1b = vec_madd(v_xb, v_mhu_A2,    v_mhu_A1);
   vdt v_sum2a = vec_madd(v_xa, v_mhu_B2,    v_mhu_B1);
   vdt v_sum2b = vec_madd(v_xb, v_mhu_B2,    v_mhu_B1);
   vdt v_sum3a = vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);
   vdt v_sum3b = vec_madd(v_xb, v_tauRdt_C2, v_tauRdt_C1);

//2
/*
   sum1[ii] = mhu_a[2]     + VM[ii]*sum1[ii];
   sum1[ii] = mhu_a[1]     + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[12]    + VM[ii]*sum2[ii];
   sum2[ii] = mhu_a[11]    + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[15] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[14] + VM[ii]*sum3[ii];
*/

   v_mhu_A1    = VEC_LDS(0, &mhu_a[1]);
   v_mhu_A2    = VEC_LDS(0, &mhu_a[2]);
   v_mhu_B1    = VEC_LDS(0, &mhu_a[11]);
   v_mhu_B2    = VEC_LDS(0, &mhu_a[12]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[14]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[15]);

   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum1b = vec_madd(v_xb, v_sum1b, v_mhu_A2);
   v_sum1b = vec_madd(v_xb, v_sum1b, v_mhu_A1);
   v_sum2b = vec_madd(v_xb, v_sum2b, v_mhu_B2);
   v_sum2b = vec_madd(v_xb, v_sum2b, v_mhu_B1); 
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C1);

//3
/*
   sum1[ii] = mhu_a[0]     + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[10]    + VM[ii]*sum2[ii];
   sum2[ii] = mhu_a[9]     + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[13] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[12] + VM[ii]*sum3[ii];
*/

   v_mhu_A1    = VEC_LDS(0, &mhu_a[0]);
   v_mhu_B1    = VEC_LDS(0, &mhu_a[9]);
   v_mhu_B2    = VEC_LDS(0, &mhu_a[10]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[12]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[13]);

   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum1b = vec_madd(v_xb, v_sum1b, v_mhu_A1);
   v_sum2b = vec_madd(v_xb, v_sum2b, v_mhu_B2);
   v_sum2b = vec_madd(v_xb, v_sum2b, v_mhu_B1); 
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C1);

//4
/*
   sum2[ii] = mhu_a[8]     + VM[ii]*sum2[ii];
   sum2[ii] = mhu_a[7]     + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[11] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[10] + VM[ii]*sum3[ii];
*/
   v_mhu_B1    = VEC_LDS(0, &mhu_a[7]);
   v_mhu_B2    = VEC_LDS(0, &mhu_a[8]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[10]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[11]);

   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum2b = vec_madd(v_xb, v_sum2b, v_mhu_B2);
   v_sum2b = vec_madd(v_xb, v_sum2b, v_mhu_B1); 
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C1);

//5
/*
   sum2[ii] =   mhu_a[6]   + VM[ii]*sum2[ii];
   sum2[ii] =   mhu_a[5]   + VM[ii]*sum2[ii];
   sum3[ii] =  tauRdt_a[9] + VM[ii]*sum3[ii];
   sum3[ii] =  tauRdt_a[8] + VM[ii]*sum3[ii];
*/
   v_mhu_B1    = VEC_LDS(0, &mhu_a[5]);
   v_mhu_B2    = VEC_LDS(0, &mhu_a[6]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[8]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[9]);

   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum2b = vec_madd(v_xb, v_sum2b, v_mhu_B2);
   v_sum2b = vec_madd(v_xb, v_sum2b, v_mhu_B1); 
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C1);

//6
/*
   sum3[ii] = tauRdt_a[7] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[6] + VM[ii]*sum3[ii];
*/
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[6]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[7]);

   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C1);

//7
/*
   sum3[ii] = tauRdt_a[5] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[4] + VM[ii]*sum3[ii];
*/
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[4]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[5]);

   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C1);

//8
/*
   sum3[ii] = tauRdt_a[3] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[2] + VM[ii]*sum3[ii];
*/
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[2]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[3]);

   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C1);

//9
/*
   sum3[ii] = tauRdt_a[1] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[0] + VM[ii]*sum3[ii];
*/
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[0]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[1]);

   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C1);



/*
   double mhu= sum1/sum2;
   for (int i=exp_l;i>=0;i--) sum4 = exp_a[i]+sum4*sum3; 
   double tauRdt =   1.0 - 1.0/sum4 ;
    g[ii] +=  mhu*tauRdt - g[ii]*tauRdt;
*/

// A
   v_mhu_A2 = vec_swdiv_nochk(v_sum1a,v_sum2a); 

       //vdt v_sum4  = vec_madd(v_sum3a, v_exp_a16,    v_exp_a15);
       //v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a14);
       //v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a13);
       //v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a12);
       //v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a11);
       vdt v_sum4  = vec_madd(v_sum3a, v_exp_a12,    v_exp_a11);
       v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a10);
       v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a9);
       v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a8);
       v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a7);
       v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a6);
       v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a5);
       v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a4);
       v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a3);
       v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a2);
       v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a1);
       v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a0);

   v_sum4 = vec_re(v_sum4);

   v_tauRdt_C2 = vec_sub(v_ONE, v_sum4);
 
   v_mhu_B2    = vec_ld(0, &g[ii]); 
 
   v_xa        = vec_sub(v_mhu_A2, v_mhu_B2); 
 
   v_mhu_B2    = vec_madd(v_tauRdt_C2, v_xa, v_mhu_B2);

   vec_st(v_mhu_B2, 0, &g[ii]);

// B

   v_mhu_A1 = vec_swdiv_nochk(v_sum1b,v_sum2b); 

   //v_sum4  = vec_madd(v_sum3b, v_exp_a16,    v_exp_a15);
   //v_sum4  = vec_madd(v_sum3b, v_sum4,    v_exp_a14);
   //v_sum4  = vec_madd(v_sum3b, v_sum4,    v_exp_a13);
   //v_sum4  = vec_madd(v_sum3b, v_sum4,    v_exp_a12);
   //v_sum4  = vec_madd(v_sum3b, v_sum4,    v_exp_a11);
   v_sum4  = vec_madd(v_sum3b, v_exp_a12,    v_exp_a11);
   v_sum4  = vec_madd(v_sum3b, v_sum4,    v_exp_a10);
   v_sum4  = vec_madd(v_sum3b, v_sum4,    v_exp_a9);
   v_sum4  = vec_madd(v_sum3b, v_sum4,    v_exp_a8);
   v_sum4  = vec_madd(v_sum3b, v_sum4,    v_exp_a7);
   v_sum4  = vec_madd(v_sum3b, v_sum4,    v_exp_a6);
   v_sum4  = vec_madd(v_sum3b, v_sum4,    v_exp_a5);
   v_sum4  = vec_madd(v_sum3b, v_sum4,    v_exp_a4);
   v_sum4  = vec_madd(v_sum3b, v_sum4,    v_exp_a3);
   v_sum4  = vec_madd(v_sum3b, v_sum4,    v_exp_a2);
   v_sum4  = vec_madd(v_sum3b, v_sum4,    v_exp_a1);
   v_sum4  = vec_madd(v_sum3b, v_sum4,    v_exp_a0);

   v_sum4 = vec_re(v_sum4);

   v_tauRdt_C1 = vec_sub(v_ONE, v_sum4);

   v_mhu_B1    = vec_ld(0, &g[ii+4]); 
 
   v_xb        = vec_sub(v_mhu_A1, v_mhu_B1); 
 
   v_mhu_B1    = vec_madd(v_tauRdt_C1, v_xb, v_mhu_B1);

   vec_st(v_mhu_B1, 0, &g[ii+4]);

  }
 for (ii;ii<TTROUND(nCells,4);ii+=4)
 { 
   //BODYSTART
   vdt v_xa = vec_ld(0, &VM[ii]);  

   vdt v_mhu_A1    = VEC_LDS(0, &mhu_a[3]);
   vdt v_mhu_A2    = VEC_LDS(0, &mhu_a[4]);
   vdt v_mhu_B1    = VEC_LDS(0, &mhu_a[13]);
   vdt v_mhu_B2    = VEC_LDS(0, &mhu_a[14]);
   vdt v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[16]);
   vdt v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[17]);

   vdt v_sum1a = vec_madd(v_xa, v_mhu_A2,    v_mhu_A1);
   vdt v_sum2a = vec_madd(v_xa, v_mhu_B2,    v_mhu_B1);
   vdt v_sum3a = vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);

//2
/*
   sum1[ii] = mhu_a[2]     + VM[ii]*sum1[ii];
   sum1[ii] = mhu_a[1]     + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[12]    + VM[ii]*sum2[ii];
   sum2[ii] = mhu_a[11]    + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[15] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[14] + VM[ii]*sum3[ii];
*/

   v_mhu_A1    = VEC_LDS(0, &mhu_a[1]);
   v_mhu_A2    = VEC_LDS(0, &mhu_a[2]);
   v_mhu_B1    = VEC_LDS(0, &mhu_a[11]);
   v_mhu_B2    = VEC_LDS(0, &mhu_a[12]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[14]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[15]);

   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);

//3
/*
   sum1[ii] = mhu_a[0]     + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[10]    + VM[ii]*sum2[ii];
   sum2[ii] = mhu_a[9]     + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[13] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[12] + VM[ii]*sum3[ii];
*/

   v_mhu_A1    = VEC_LDS(0, &mhu_a[0]);
   v_mhu_B1    = VEC_LDS(0, &mhu_a[9]);
   v_mhu_B2    = VEC_LDS(0, &mhu_a[10]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[12]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[13]);

   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);

//4
/*
   sum2[ii] = mhu_a[8]     + VM[ii]*sum2[ii];
   sum2[ii] = mhu_a[7]     + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[11] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[10] + VM[ii]*sum3[ii];
*/
   v_mhu_B1    = VEC_LDS(0, &mhu_a[7]);
   v_mhu_B2    = VEC_LDS(0, &mhu_a[8]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[10]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[11]);

   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);

//5
/*
   sum2[ii] =   mhu_a[6]   + VM[ii]*sum2[ii];
   sum2[ii] =   mhu_a[5]   + VM[ii]*sum2[ii];
   sum3[ii] =  tauRdt_a[9] + VM[ii]*sum3[ii];
   sum3[ii] =  tauRdt_a[8] + VM[ii]*sum3[ii];
*/
   v_mhu_B1    = VEC_LDS(0, &mhu_a[5]);
   v_mhu_B2    = VEC_LDS(0, &mhu_a[6]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[8]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[9]);

   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);

//6
/*
   sum3[ii] = tauRdt_a[7] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[6] + VM[ii]*sum3[ii];
*/
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[6]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[7]);

   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);

//7
/*
   sum3[ii] = tauRdt_a[5] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[4] + VM[ii]*sum3[ii];
*/
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[4]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[5]);

   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);

//8
/*
   sum3[ii] = tauRdt_a[3] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[2] + VM[ii]*sum3[ii];
*/
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[2]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[3]);

   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);

//9
/*
   sum3[ii] = tauRdt_a[1] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[0] + VM[ii]*sum3[ii];
*/
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[0]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[1]);

   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);



/*
   double mhu= sum1/sum2;
   for (int i=exp_l;i>=0;i--) sum4 = exp_a[i]+sum4*sum3; 
   double tauRdt =   1.0 - 1.0/sum4 ;
    g[ii] +=  mhu*tauRdt - g[ii]*tauRdt;
*/

// A
   v_mhu_A2 = vec_swdiv_nochk(v_sum1a,v_sum2a); 

       //vdt v_sum4  = vec_madd(v_sum3a, v_exp_a16,    v_exp_a15);
       //v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a14);
       //v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a13);
       //v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a12);
       //v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a11);
       vdt v_sum4  = vec_madd(v_sum3a, v_exp_a12,    v_exp_a11);
       v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a10);
       v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a9);
       v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a8);
       v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a7);
       v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a6);
       v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a5);
       v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a4);
       v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a3);
       v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a2);
       v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a1);
       v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a0);

   v_sum4 = vec_re(v_sum4);

   v_tauRdt_C2 = vec_sub(v_ONE, v_sum4);
 
   v_mhu_B2    = vec_ld(0, &g[ii]); 
 
   v_xa        = vec_sub(v_mhu_A2, v_mhu_B2); 
 
   v_mhu_B2    = vec_madd(v_tauRdt_C2, v_xa, v_mhu_B2);

   vec_st(v_mhu_B2, 0, &g[ii]); //BODYEND

  }



}





void update_hGate_v1(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a)
{
  //update_hGate(dt, nCells, VM, g, mhu_a, tauR_a); return;
  typedef vector4double vdt;

int gateIndex=1;
int mhu_l=7 ;
int mhu_m= 9;
int tauR_l=11;
int tauR_m=11;
   
 double  tauRdt_a[11];

         tauRdt_a[0]  = tauR_a[0]*dt;
         tauRdt_a[1]  = tauR_a[1]*dt;
         tauRdt_a[2]  = tauR_a[2]*dt;
         tauRdt_a[3]  = tauR_a[3]*dt;
         tauRdt_a[4]  = tauR_a[4]*dt;
         tauRdt_a[5]  = tauR_a[5]*dt;
         tauRdt_a[6]  = tauR_a[6]*dt; 
         tauRdt_a[7]  = tauR_a[7]*dt;
         tauRdt_a[8]  = tauR_a[8]*dt;
         tauRdt_a[9]  = tauR_a[9]*dt;
         tauRdt_a[10] = tauR_a[10]*dt;

 int ii=0;

 // 0 aligned
 // 1: load produces (JUNK, &VM[0], &VM[1], &VM[2])
 // 2: load produces (JUNK, JUNK, &VM[0], &VM[1])
 // 3: load produces (JUNK, JUNK, JUNK, &VM[0])
 int misalignment = ((long int)&VM[0] % 32) / 8;
 assert(misalignment == 0);

#if 1
 if (misalignment) { 
      vdt v_xa   = vec_ld(0, &VM[ii]);  //BODYCOM0 
//BODYCOM0 
// 1//BODYCOM0 
  vdt v_mhu_A1  = VEC_LDS(0, &mhu_a[7]);//BODYCOM0 
  vdt v_mhu_A2  = VEC_LDS(0, &mhu_a[8]);//BODYCOM0 
//BODYCOM0 
  vdt v_mhu_B1 = VEC_LDS(0, &mhu_a[14]);//BODYCOM0 
  vdt v_mhu_B2 = VEC_LDS(0, &mhu_a[15]);//BODYCOM0 
  vdt v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[9]);//BODYCOM0 
  vdt v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[10]);//BODYCOM0 
//BODYCOM0 
  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[20]);//BODYCOM0 
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[21]);//BODYCOM0 
   vdt v_sum1a =  vec_madd(v_xa, v_mhu_A2, v_mhu_A1);//BODYCOM0 
   vdt v_sum2a =  vec_madd(v_xa, v_mhu_B2, v_mhu_B1);//BODYCOM0 
   vdt v_sum3a =  vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);//BODYCOM0 
   vdt v_sum4a =  vec_madd(v_xa, v_tauR_D2, v_tauR_D1);//BODYCOM0 
// 2//BODYCOM0 
   v_mhu_A1  = VEC_LDS(0, &mhu_a[5]);//BODYCOM0 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[6]);//BODYCOM0 
   v_mhu_B1 = VEC_LDS(0, &mhu_a[12]);//BODYCOM0 
   v_mhu_B2 = VEC_LDS(0, &mhu_a[13]);//BODYCOM0 
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[7]);//BODYCOM0 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[8]);//BODYCOM0 
  v_tauR_D1   = VEC_LDS(0, &tauR_a[18]);//BODYCOM0 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[19]);//BODYCOM0 
//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM0 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM0 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM0 
// 3//BODYCOM0 
//BODYCOM0 
   v_mhu_A1  = VEC_LDS(0, &mhu_a[3]);//BODYCOM0 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[4]);//BODYCOM0 
   v_mhu_B1 = VEC_LDS(0, &mhu_a[10]);//BODYCOM0 
   v_mhu_B2 = VEC_LDS(0, &mhu_a[11]);//BODYCOM0 
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[5]);//BODYCOM0 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[6]);//BODYCOM0 
  v_tauR_D1   = VEC_LDS(0, &tauR_a[16]);//BODYCOM0 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[17]);//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM0 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM0 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM0 
// 4//BODYCOM0 
//BODYCOM0 
   v_mhu_A1  = VEC_LDS(0, &mhu_a[1]);//BODYCOM0 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[2]);//BODYCOM0 
   v_mhu_B2  = VEC_LDS(0, &mhu_a[9]);//BODYCOM0 
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[3]);//BODYCOM0 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[4]);//BODYCOM0 
  v_tauR_D1   = VEC_LDS(0, &tauR_a[14]);//BODYCOM0 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[15]);//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM0 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM0 
   // v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM0 
// 5//BODYCOM0 
//BODYCOM0 
  v_mhu_A2  = VEC_LDS(0, &mhu_a[0]);//BODYCOM0 
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[1]);//BODYCOM0 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[2]);//BODYCOM0 
  v_tauR_D1   = VEC_LDS(0, &tauR_a[12]);//BODYCOM0 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[13]);//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM0 
   // v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM0 
   // v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM0 
   // v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM0 
// 6//BODYCOM0 
//BODYCOM0 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[0]);//BODYCOM0 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[11]);//BODYCOM0 
   // v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM0 
   // v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM0 
   // v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM0 
   // v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   //v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM0 
   //v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM0 
   // v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A2);//BODYCOM0 
   // v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A1);//BODYCOM0 
   // v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B2);//BODYCOM0 
   // v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B1); //BODYCOM0 
//BODYCOM0 
// Still need to reuse the register I have above ... the code will look like gibberish because//BODYCOM0 
// the targets won't make sense.//BODYCOM0 
// A//BODYCOM0 
   v_xa = vec_swdiv_nochk(v_sum3a,v_sum4a); //BODYCOM0 
   v_mhu_A2    = vec_swdiv_nochk(v_sum1a,v_sum2a); //BODYCOM0 
   v_mhu_B2   = vec_ld(0, &g[ii]);  //BODYCOM0 
   v_tauR_D2 = vec_sub(v_mhu_A2, v_mhu_B2);  //BODYCOM0 
   v_mhu_B2       = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);//BODYCOM0 

   
   //vec_st(v_mhu_B2, 0, &g[ii]);
   switch (misalignment) {
   case 1:
     g[ii+0] = v_mhu_B2 get [1];
     g[ii+1] = v_mhu_B2 get [2];
     g[ii+2] = v_mhu_B2 get [3];
     break;

   case 2:
     g[ii+0] = v_mhu_B2 get [2];
     g[ii+1] = v_mhu_B2 get [3];
     break;

   case 3:
     g[ii+0] = v_mhu_B2 get [3];
     break;

     default:
    //println("BAD MISALIGNMENT: %d\n", misalignment); fflush(stdout);
     assert(0);
   }

   ii = 4 - misalignment;
 }

 for (ii;ii<TTROUND(nCells,8);ii+=8)
 { 

   vdt v_xa   = vec_ld(0, &VM[ii]);  
   vdt v_xb   = vec_ld(0, &VM[ii+4]);  

// 1
  vdt v_mhu_A1  = VEC_LDS(0, &mhu_a[7]);
  vdt v_mhu_A2  = VEC_LDS(0, &mhu_a[8]);

  vdt v_mhu_B1 = VEC_LDS(0, &mhu_a[14]);
  vdt v_mhu_B2 = VEC_LDS(0, &mhu_a[15]);
  vdt v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[9]);
  vdt v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[10]);

  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[20]);
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[21]);
   vdt v_sum1a =  vec_madd(v_xa, v_mhu_A2, v_mhu_A1);
   vdt v_sum1b =  vec_madd(v_xb, v_mhu_A2, v_mhu_A1);

   vdt v_sum2a =  vec_madd(v_xa, v_mhu_B2, v_mhu_B1);
   vdt v_sum2b =  vec_madd(v_xb, v_mhu_B2, v_mhu_B1);

   vdt v_sum3a =  vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);
   vdt v_sum3b =  vec_madd(v_xb, v_tauRdt_C2, v_tauRdt_C1);

   vdt v_sum4a =  vec_madd(v_xa, v_tauR_D2, v_tauR_D1);
   vdt v_sum4b =  vec_madd(v_xb, v_tauR_D2, v_tauR_D1);
// 2
   v_mhu_A1  = VEC_LDS(0, &mhu_a[5]);
   v_mhu_A2  = VEC_LDS(0, &mhu_a[6]);
   v_mhu_B1 = VEC_LDS(0, &mhu_a[12]);
   v_mhu_B2 = VEC_LDS(0, &mhu_a[13]);
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[7]);
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[8]);
  v_tauR_D1   = VEC_LDS(0, &tauR_a[18]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[19]);

   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A2);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A1);
   v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B2);
   v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B1); 
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D2);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D1);
// 3

   v_mhu_A1  = VEC_LDS(0, &mhu_a[3]);
   v_mhu_A2  = VEC_LDS(0, &mhu_a[4]);
   v_mhu_B1 = VEC_LDS(0, &mhu_a[10]);
   v_mhu_B2 = VEC_LDS(0, &mhu_a[11]);
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[5]);
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[6]);
  v_tauR_D1   = VEC_LDS(0, &tauR_a[16]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[17]);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A2);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A1);
   v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B2);
   v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B1); 
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D2);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D1);
// 4

   v_mhu_A1  = VEC_LDS(0, &mhu_a[1]);
   v_mhu_A2  = VEC_LDS(0, &mhu_a[2]);
   v_mhu_B2  = VEC_LDS(0, &mhu_a[9]);
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[3]);
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[4]);
  v_tauR_D1   = VEC_LDS(0, &tauR_a[14]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[15]);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   // v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A2);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A1);
   v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B2);
   // v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B1); 
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D2);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D1);
// 5

  v_mhu_A2  = VEC_LDS(0, &mhu_a[0]);
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[1]);
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[2]);
  v_tauR_D1   = VEC_LDS(0, &tauR_a[12]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[13]);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   // v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   // v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   // v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A2);
   // v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A1);
   // v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B2);
   // v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B1); 
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D2);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D1);
// 6

  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[0]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[11]);
   // v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   // v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   // v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   // v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   //v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   //v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);
   // v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A2);
   // v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A1);
   // v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B2);
   // v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B1); 
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   //v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D2);
   //v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D1);

// Still need to reuse the register I have above ... the code will look like gibberish because
// the targets won't make sense.
// A
   v_xa = vec_swdiv_nochk(v_sum3a,v_sum4a); 
   v_mhu_A2    = vec_swdiv_nochk(v_sum1a,v_sum2a); 
   v_mhu_B2   = vec_ld(0, &g[ii]);  
   v_tauR_D2 = vec_sub(v_mhu_A2, v_mhu_B2);  
   v_mhu_B2       = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);
   vec_st(v_mhu_B2, 0, &g[ii]);
// B
   v_xb = vec_swdiv_nochk(v_sum3b,v_sum4b); 
   v_mhu_A1    = vec_swdiv_nochk(v_sum1b,v_sum2b); 
   v_mhu_B1   = vec_ld(0, &g[ii+4]);  
   v_tauR_D2 = vec_sub(v_mhu_A1, v_mhu_B1);  
   v_mhu_B1       = vec_madd(v_tauR_D2, v_xb, v_mhu_B1);
   vec_st(v_mhu_B1, 0, &g[ii+4]);
 }
#endif

 for (ii;ii<TTROUND(nCells,4);ii+=4)
 { 
   //BODYSTART
   vdt v_xa   = vec_ld(0, &VM[ii]);  

// 1
  vdt v_mhu_A1  = VEC_LDS(0, &mhu_a[7]);
  vdt v_mhu_A2  = VEC_LDS(0, &mhu_a[8]);

  vdt v_mhu_B1 = VEC_LDS(0, &mhu_a[14]);
  vdt v_mhu_B2 = VEC_LDS(0, &mhu_a[15]);
  vdt v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[9]);
  vdt v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[10]);

  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[20]);
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[21]);
   vdt v_sum1a =  vec_madd(v_xa, v_mhu_A2, v_mhu_A1);
   vdt v_sum2a =  vec_madd(v_xa, v_mhu_B2, v_mhu_B1);
   vdt v_sum3a =  vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);
   vdt v_sum4a =  vec_madd(v_xa, v_tauR_D2, v_tauR_D1);
// 2
   v_mhu_A1  = VEC_LDS(0, &mhu_a[5]);
   v_mhu_A2  = VEC_LDS(0, &mhu_a[6]);
   v_mhu_B1 = VEC_LDS(0, &mhu_a[12]);
   v_mhu_B2 = VEC_LDS(0, &mhu_a[13]);
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[7]);
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[8]);
  v_tauR_D1   = VEC_LDS(0, &tauR_a[18]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[19]);

   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);
// 3

   v_mhu_A1  = VEC_LDS(0, &mhu_a[3]);
   v_mhu_A2  = VEC_LDS(0, &mhu_a[4]);
   v_mhu_B1 = VEC_LDS(0, &mhu_a[10]);
   v_mhu_B2 = VEC_LDS(0, &mhu_a[11]);
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[5]);
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[6]);
  v_tauR_D1   = VEC_LDS(0, &tauR_a[16]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[17]);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);
// 4

   v_mhu_A1  = VEC_LDS(0, &mhu_a[1]);
   v_mhu_A2  = VEC_LDS(0, &mhu_a[2]);
   v_mhu_B2  = VEC_LDS(0, &mhu_a[9]);
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[3]);
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[4]);
  v_tauR_D1   = VEC_LDS(0, &tauR_a[14]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[15]);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   // v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);
// 5

  v_mhu_A2  = VEC_LDS(0, &mhu_a[0]);
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[1]);
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[2]);
  v_tauR_D1   = VEC_LDS(0, &tauR_a[12]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[13]);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   // v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   // v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   // v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);
// 6

  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[0]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[11]);
   // v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   // v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   // v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   // v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   //v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   //v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);
   // v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A2);
   // v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A1);
   // v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B2);
   // v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B1); 

// Still need to reuse the register I have above ... the code will look like gibberish because
// the targets won't make sense.
// A
   v_xa = vec_swdiv_nochk(v_sum3a,v_sum4a); 
   v_mhu_A2    = vec_swdiv_nochk(v_sum1a,v_sum2a); 
   v_mhu_B2   = vec_ld(0, &g[ii]);  
   v_tauR_D2 = vec_sub(v_mhu_A2, v_mhu_B2);  
   v_mhu_B2       = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);
   vec_st(v_mhu_B2, 0, &g[ii]); //BODYEND
 }

#if 1
 // the elements in the end of the vector should be ignored
 int leftover = nCells - ii;
 if (leftover) { 
      vdt v_xa   = vec_ld(0, &VM[ii]);  //BODYCOM1 
//BODYCOM1 
// 1//BODYCOM1 
  vdt v_mhu_A1  = VEC_LDS(0, &mhu_a[7]);//BODYCOM1 
  vdt v_mhu_A2  = VEC_LDS(0, &mhu_a[8]);//BODYCOM1 
//BODYCOM1 
  vdt v_mhu_B1 = VEC_LDS(0, &mhu_a[14]);//BODYCOM1 
  vdt v_mhu_B2 = VEC_LDS(0, &mhu_a[15]);//BODYCOM1 
  vdt v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[9]);//BODYCOM1 
  vdt v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[10]);//BODYCOM1 
//BODYCOM1 
  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[20]);//BODYCOM1 
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[21]);//BODYCOM1 
   vdt v_sum1a =  vec_madd(v_xa, v_mhu_A2, v_mhu_A1);//BODYCOM1 
   vdt v_sum2a =  vec_madd(v_xa, v_mhu_B2, v_mhu_B1);//BODYCOM1 
   vdt v_sum3a =  vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);//BODYCOM1 
   vdt v_sum4a =  vec_madd(v_xa, v_tauR_D2, v_tauR_D1);//BODYCOM1 
// 2//BODYCOM1 
   v_mhu_A1  = VEC_LDS(0, &mhu_a[5]);//BODYCOM1 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[6]);//BODYCOM1 
   v_mhu_B1 = VEC_LDS(0, &mhu_a[12]);//BODYCOM1 
   v_mhu_B2 = VEC_LDS(0, &mhu_a[13]);//BODYCOM1 
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[7]);//BODYCOM1 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[8]);//BODYCOM1 
  v_tauR_D1   = VEC_LDS(0, &tauR_a[18]);//BODYCOM1 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[19]);//BODYCOM1 
//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM1 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM1 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM1 
// 3//BODYCOM1 
//BODYCOM1 
   v_mhu_A1  = VEC_LDS(0, &mhu_a[3]);//BODYCOM1 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[4]);//BODYCOM1 
   v_mhu_B1 = VEC_LDS(0, &mhu_a[10]);//BODYCOM1 
   v_mhu_B2 = VEC_LDS(0, &mhu_a[11]);//BODYCOM1 
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[5]);//BODYCOM1 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[6]);//BODYCOM1 
  v_tauR_D1   = VEC_LDS(0, &tauR_a[16]);//BODYCOM1 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[17]);//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM1 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM1 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM1 
// 4//BODYCOM1 
//BODYCOM1 
   v_mhu_A1  = VEC_LDS(0, &mhu_a[1]);//BODYCOM1 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[2]);//BODYCOM1 
   v_mhu_B2  = VEC_LDS(0, &mhu_a[9]);//BODYCOM1 
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[3]);//BODYCOM1 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[4]);//BODYCOM1 
  v_tauR_D1   = VEC_LDS(0, &tauR_a[14]);//BODYCOM1 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[15]);//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM1 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM1 
   // v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM1 
// 5//BODYCOM1 
//BODYCOM1 
  v_mhu_A2  = VEC_LDS(0, &mhu_a[0]);//BODYCOM1 
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[1]);//BODYCOM1 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[2]);//BODYCOM1 
  v_tauR_D1   = VEC_LDS(0, &tauR_a[12]);//BODYCOM1 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[13]);//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM1 
   // v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM1 
   // v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM1 
   // v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM1 
// 6//BODYCOM1 
//BODYCOM1 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[0]);//BODYCOM1 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[11]);//BODYCOM1 
   // v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM1 
   // v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM1 
   // v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM1 
   // v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   //v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM1 
   //v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM1 
   // v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A2);//BODYCOM1 
   // v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A1);//BODYCOM1 
   // v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B2);//BODYCOM1 
   // v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B1); //BODYCOM1 
//BODYCOM1 
// Still need to reuse the register I have above ... the code will look like gibberish because//BODYCOM1 
// the targets won't make sense.//BODYCOM1 
// A//BODYCOM1 
   v_xa = vec_swdiv_nochk(v_sum3a,v_sum4a); //BODYCOM1 
   v_mhu_A2    = vec_swdiv_nochk(v_sum1a,v_sum2a); //BODYCOM1 
   v_mhu_B2   = vec_ld(0, &g[ii]);  //BODYCOM1 
   v_tauR_D2 = vec_sub(v_mhu_A2, v_mhu_B2);  //BODYCOM1 
   v_mhu_B2       = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);//BODYCOM1 

   
   //vec_st(v_mhu_B2, 0, &g[ii]);
   switch (leftover) {
   case 1:
     g[ii+0] = v_mhu_B2 get [0];
     break;

   case 2:
     g[ii+0] = v_mhu_B2 get [0];
     g[ii+1] = v_mhu_B2 get [1];
     break;

   case 3:
     g[ii+0] = v_mhu_B2 get [0];
     g[ii+1] = v_mhu_B2 get [1];
     g[ii+2] = v_mhu_B2 get [2];
     break;

     default:
    //println("BAD LEFTOVER: %d\n", leftover); fflush(stdout);
    assert(0);
   }
 }
#endif
}





void update_jGate_v1(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a)
{
  //update_jGate(dt, nCells, VM, g, mhu_a, tauR_a); return;
  typedef vector4double vdt;

int gateIndex=2;
int mhu_l=7 ;
int mhu_m= 9;
int tauR_l=1;
int tauR_m=13;
   
 double  tauRdt_a[13];


         tauRdt_a[0]  = tauR_a[0]*dt;
         tauRdt_a[1]  = tauR_a[1]*dt;
         tauRdt_a[2]  = tauR_a[2]*dt;
         tauRdt_a[3]  = tauR_a[3]*dt;
         tauRdt_a[4]  = tauR_a[4]*dt;
         tauRdt_a[5]  = tauR_a[5]*dt;
         tauRdt_a[6]  = tauR_a[6]*dt; 
         tauRdt_a[7]  = tauR_a[7]*dt;
         tauRdt_a[8]  = tauR_a[8]*dt;
         tauRdt_a[9]  = tauR_a[9]*dt;
         tauRdt_a[10] = tauR_a[10]*dt;
         tauRdt_a[11] = tauR_a[11]*dt;
         tauRdt_a[12] = tauR_a[12]*dt;
 int ii=0;

 // 0 aligned
 // 1: load produces (JUNK, &VM[0], &VM[1], &VM[2])
 // 2: load produces (JUNK, JUNK, &VM[0], &VM[1])
 // 3: load produces (JUNK, JUNK, JUNK, &VM[0])
 int misalignment = ((long int)&VM[0] % 32) / 8;
 assert(misalignment == 0);

 if (misalignment) { 
      vdt v_xa   = vec_ld(0, &VM[ii]);  //BODYCOM0 
//BODYCOM0 
// 1//BODYCOM0 
   vdt v_mhu_A1  = VEC_LDS(0, &mhu_a[7]);//BODYCOM0 
   vdt v_mhu_A2  = VEC_LDS(0, &mhu_a[8]);//BODYCOM0 
   vdt v_mhu_B1 = VEC_LDS(0, &mhu_a[14]);//BODYCOM0 
   vdt v_mhu_B2 = VEC_LDS(0, &mhu_a[15]);//BODYCOM0 
   vdt v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[11]);//BODYCOM0 
   vdt v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[12]);//BODYCOM0 
   vdt v_sum1a =  vec_madd(v_xa, v_mhu_A2, v_mhu_A1);//BODYCOM0 
   vdt v_sum2a =  vec_madd(v_xa, v_mhu_B2, v_mhu_B1);//BODYCOM0 
   vdt v_sum3a =  vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);//BODYCOM0 
// 2//BODYCOM0 
   v_mhu_A1  = VEC_LDS(0, &mhu_a[5]);//BODYCOM0 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[6]);//BODYCOM0 
   v_mhu_B1 = VEC_LDS(0, &mhu_a[12]);//BODYCOM0 
   v_mhu_B2 = VEC_LDS(0, &mhu_a[13]);//BODYCOM0 
   v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[9]);//BODYCOM0 
   v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[10]);//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM0 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM0 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
// 3//BODYCOM0 
   v_mhu_A1  = VEC_LDS(0, &mhu_a[3]);//BODYCOM0 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[4]);//BODYCOM0 
   v_mhu_B1 = VEC_LDS(0, &mhu_a[10]);//BODYCOM0 
   v_mhu_B2 = VEC_LDS(0, &mhu_a[11]);//BODYCOM0 
   v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[7]);//BODYCOM0 
   v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[8]);//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM0 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM0 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
// 4//BODYCOM0 
   v_mhu_A1  = VEC_LDS(0, &mhu_a[1]);//BODYCOM0 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[2]);//BODYCOM0 
   v_mhu_B2  = VEC_LDS(0, &mhu_a[9]);//BODYCOM0 
   v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[5]);//BODYCOM0 
   v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[6]);//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM0 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
// 5//BODYCOM0 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[0]);//BODYCOM0 
   v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[3]);//BODYCOM0 
   v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[4]);//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
// 6//BODYCOM0 
   v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[1]);//BODYCOM0 
   v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[2]);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
// 7//BODYCOM0 
   v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[0]);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
//BODYCOM0 
// 1st//BODYCOM0 
   vdt temp1a    = vec_swdiv_nochk(v_sum1a,v_sum2a); //BODYCOM0 
   vdt temp2a   = vec_ld(0, &g[ii]);  //BODYCOM0 
   temp1a = vec_sub(temp1a, temp2a);  //BODYCOM0 
   temp2a       = vec_madd(temp1a, v_sum3a, temp2a);//BODYCOM0 

   
   //vec_st(temp2a, 0, &g[ii]);
   switch (misalignment) {
   case 1:
     g[ii+0] = temp2a get [1];
     g[ii+1] = temp2a get [2];
     g[ii+2] = temp2a get [3];
     break;

   case 2:
     g[ii+0] = temp2a get [2];
     g[ii+1] = temp2a get [3];
     break;

   case 3:
     g[ii+0] = temp2a get [3];
     break;

     default:
    //println("BAD MISALIGNMENT: %d\n", misalignment); fflush(stdout);
     assert(0);
   }

   ii = 4 - misalignment;
 }

 //for (ii=0;ii<0;ii+=16)
 for (ii;ii<TTROUND(nCells,16);ii+=16)
 { 
   vdt v_xa   = vec_ld(0, &VM[ii]);  
   vdt v_xb   = vec_ld(0, &VM[ii+4]);  
   vdt v_xc   = vec_ld(0, &VM[ii+8]);  
   vdt v_xd   = vec_ld(0, &VM[ii+12]);  

// 1
   vdt v_mhu_A1  = VEC_LDS(0, &mhu_a[7]);
   vdt v_mhu_A2  = VEC_LDS(0, &mhu_a[8]);
   vdt v_mhu_B1 = VEC_LDS(0, &mhu_a[14]);
   vdt v_mhu_B2 = VEC_LDS(0, &mhu_a[15]);
   vdt v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[11]);
   vdt v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[12]);
   vdt v_sum1a =  vec_madd(v_xa, v_mhu_A2, v_mhu_A1);
   vdt v_sum1b =  vec_madd(v_xb, v_mhu_A2, v_mhu_A1);
   vdt v_sum1c =  vec_madd(v_xc, v_mhu_A2, v_mhu_A1);
   vdt v_sum1d =  vec_madd(v_xd, v_mhu_A2, v_mhu_A1);
   vdt v_sum2a =  vec_madd(v_xa, v_mhu_B2, v_mhu_B1);
   vdt v_sum2b =  vec_madd(v_xb, v_mhu_B2, v_mhu_B1);
   vdt v_sum2c =  vec_madd(v_xc, v_mhu_B2, v_mhu_B1);
   vdt v_sum2d =  vec_madd(v_xd, v_mhu_B2, v_mhu_B1);
   vdt v_sum3a =  vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);
   vdt v_sum3b =  vec_madd(v_xb, v_tauRdt_C2, v_tauRdt_C1);
   vdt v_sum3c =  vec_madd(v_xc, v_tauRdt_C2, v_tauRdt_C1);
   vdt v_sum3d =  vec_madd(v_xd, v_tauRdt_C2, v_tauRdt_C1);
// 2
   v_mhu_A1  = VEC_LDS(0, &mhu_a[5]);
   v_mhu_A2  = VEC_LDS(0, &mhu_a[6]);
   v_mhu_B1 = VEC_LDS(0, &mhu_a[12]);
   v_mhu_B2 = VEC_LDS(0, &mhu_a[13]);
   v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[9]);
   v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[10]);

   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A2);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A1);
   v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B2);
   v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B1); 
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   v_sum1c =  vec_madd(v_xc, v_sum1c, v_mhu_A2);
   v_sum1c =  vec_madd(v_xc, v_sum1c, v_mhu_A1);
   v_sum2c =  vec_madd(v_xc, v_sum2c, v_mhu_B2);
   v_sum2c =  vec_madd(v_xc, v_sum2c, v_mhu_B1); 
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C2);
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C1);
   v_sum1d =  vec_madd(v_xd, v_sum1d, v_mhu_A2);
   v_sum1d =  vec_madd(v_xd, v_sum1d, v_mhu_A1);
   v_sum2d =  vec_madd(v_xd, v_sum2d, v_mhu_B2);
   v_sum2d =  vec_madd(v_xd, v_sum2d, v_mhu_B1); 
   v_sum3d =  vec_madd(v_xd, v_sum3d, v_tauRdt_C2);
   v_sum3d =  vec_madd(v_xd, v_sum3d, v_tauRdt_C1);
// 3
   v_mhu_A1  = VEC_LDS(0, &mhu_a[3]);
   v_mhu_A2  = VEC_LDS(0, &mhu_a[4]);
   v_mhu_B1 = VEC_LDS(0, &mhu_a[10]);
   v_mhu_B2 = VEC_LDS(0, &mhu_a[11]);
   v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[7]);
   v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[8]);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A2);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A1);
   v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B2);
   v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B1); 
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   v_sum1c =  vec_madd(v_xc, v_sum1c, v_mhu_A2);
   v_sum1c =  vec_madd(v_xc, v_sum1c, v_mhu_A1);
   v_sum2c =  vec_madd(v_xc, v_sum2c, v_mhu_B2);
   v_sum2c =  vec_madd(v_xc, v_sum2c, v_mhu_B1); 
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C2);
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C1);
   v_sum1d =  vec_madd(v_xd, v_sum1d, v_mhu_A2);
   v_sum1d =  vec_madd(v_xd, v_sum1d, v_mhu_A1);
   v_sum2d =  vec_madd(v_xd, v_sum2d, v_mhu_B2);
   v_sum2d =  vec_madd(v_xd, v_sum2d, v_mhu_B1); 
   v_sum3d =  vec_madd(v_xd, v_sum3d, v_tauRdt_C2);
   v_sum3d =  vec_madd(v_xd, v_sum3d, v_tauRdt_C1);
// 4
   v_mhu_A1  = VEC_LDS(0, &mhu_a[1]);
   v_mhu_A2  = VEC_LDS(0, &mhu_a[2]);
   v_mhu_B2  = VEC_LDS(0, &mhu_a[9]);
   v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[5]);
   v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[6]);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   // v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A2);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A1);
   v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B2);
   // v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B1); 
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   v_sum1c =  vec_madd(v_xc, v_sum1c, v_mhu_A2);
   v_sum1c =  vec_madd(v_xc, v_sum1c, v_mhu_A1);
   v_sum2c =  vec_madd(v_xc, v_sum2c, v_mhu_B2);
   // v_sum2c =  vec_madd(v_xc, v_sum2c, v_mhu_B1); 
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C2);
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C1);
   v_sum1d =  vec_madd(v_xd, v_sum1d, v_mhu_A2);
   v_sum1d =  vec_madd(v_xd, v_sum1d, v_mhu_A1);
   v_sum2d =  vec_madd(v_xd, v_sum2d, v_mhu_B2);
   // v_sum2d =  vec_madd(v_xd, v_sum2d, v_mhu_B1); 
   v_sum3d =  vec_madd(v_xd, v_sum3d, v_tauRdt_C2);
   v_sum3d =  vec_madd(v_xd, v_sum3d, v_tauRdt_C1);
// 5
   v_mhu_A2  = VEC_LDS(0, &mhu_a[0]);
   v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[3]);
   v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[4]);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   // v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   // v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   // v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A2);
   // v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A1);
   // v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B2);
   // v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B1); 
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   v_sum1c =  vec_madd(v_xc, v_sum1c, v_mhu_A2);
   // v_sum1c =  vec_madd(v_xc, v_sum1c, v_mhu_A1);
   // v_sum2c =  vec_madd(v_xc, v_sum2c, v_mhu_B2);
   // v_sum2c =  vec_madd(v_xc, v_sum2c, v_mhu_B1); 
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C2);
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C1);
   v_sum1d =  vec_madd(v_xd, v_sum1d, v_mhu_A2);
   // v_sum1d =  vec_madd(v_xd, v_sum1d, v_mhu_A1);
   // v_sum2d =  vec_madd(v_xd, v_sum2d, v_mhu_B2);
   // v_sum2d =  vec_madd(v_xd, v_sum2d, v_mhu_B1); 
   v_sum3d =  vec_madd(v_xd, v_sum3d, v_tauRdt_C2);
   v_sum3d =  vec_madd(v_xd, v_sum3d, v_tauRdt_C1);
// 6
   v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[1]);
   v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[2]);
   // v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   // v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   // v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   // v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   // v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A2);
   // v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A1);
   // v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B2);
   // v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B1); 
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   // v_sum1c =  vec_madd(v_xc, v_sum1c, v_mhu_A2);
   // v_sum1c =  vec_madd(v_xc, v_sum1c, v_mhu_A1);
   // v_sum2c =  vec_madd(v_xc, v_sum2c, v_mhu_B2);
   // v_sum2c =  vec_madd(v_xc, v_sum2c, v_mhu_B1); 
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C2);
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C1);
   // v_sum1d =  vec_madd(v_xd, v_sum1d, v_mhu_A2);
   // v_sum1d =  vec_madd(v_xd, v_sum1d, v_mhu_A1);
   // v_sum2d =  vec_madd(v_xd, v_sum2d, v_mhu_B2);
   // v_sum2d =  vec_madd(v_xd, v_sum2d, v_mhu_B1); 
   v_sum3d =  vec_madd(v_xd, v_sum3d, v_tauRdt_C2);
   v_sum3d =  vec_madd(v_xd, v_sum3d, v_tauRdt_C1);
// 7
   v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[0]);
   // v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   // v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   // v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   // v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   //v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   // v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A2);
   // v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A1);
   // v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B2);
   // v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B1); 
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   //v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   // v_sum1c =  vec_madd(v_xc, v_sum1c, v_mhu_A2);
   // v_sum1c =  vec_madd(v_xc, v_sum1c, v_mhu_A1);
   // v_sum2c =  vec_madd(v_xc, v_sum2c, v_mhu_B2);
   // v_sum2c =  vec_madd(v_xc, v_sum2c, v_mhu_B1); 
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C2);
   //v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C1);
   // v_sum1d =  vec_madd(v_xd, v_sum1d, v_mhu_A2);
   // v_sum1d =  vec_madd(v_xd, v_sum1d, v_mhu_A1);
   // v_sum2d =  vec_madd(v_xd, v_sum2d, v_mhu_B2);
   // v_sum2d =  vec_madd(v_xd, v_sum2d, v_mhu_B1); 
   v_sum3d =  vec_madd(v_xd, v_sum3d, v_tauRdt_C2);
   //v_sum3d =  vec_madd(v_xd, v_sum3d, v_tauRdt_C1);

// 1st
   vdt temp1a    = vec_swdiv_nochk(v_sum1a,v_sum2a); 
   vdt temp2a   = vec_ld(0, &g[ii]);  
   temp1a = vec_sub(temp1a, temp2a);  
   temp2a       = vec_madd(temp1a, v_sum3a, temp2a);
   vec_st(temp2a, 0, &g[ii]);
// 2nd
   vdt temp1b    = vec_swdiv_nochk(v_sum1b,v_sum2b); 
   vdt temp2b   = vec_ld(0, &g[ii+4]);  
   temp1b = vec_sub(temp1b, temp2b);  
   temp2b       = vec_madd(temp1b, v_sum3b, temp2b);
   vec_st(temp2b, 0, &g[ii+4]);
// 3rd
   vdt temp1c    = vec_swdiv_nochk(v_sum1c,v_sum2c); 
   vdt temp2c   = vec_ld(0, &g[ii+8]);  
   temp1c = vec_sub(temp1c, temp2c);  
   temp2c       = vec_madd(temp1c, v_sum3c, temp2c);
   vec_st(temp2c, 0, &g[ii+8]);
// 4th 
   vdt temp1d    = vec_swdiv_nochk(v_sum1d,v_sum2d); 
   vdt temp2d   = vec_ld(0, &g[ii+12]);  
   temp1d = vec_sub(temp1d, temp2d);  
   temp2d       = vec_madd(temp1d, v_sum3d, temp2d);
   vec_st(temp2d, 0, &g[ii+12]);
 }
 for (ii;ii<TTROUND(nCells,4);ii+=4)
 { 
   //BODYSTART
   vdt v_xa   = vec_ld(0, &VM[ii]);  

// 1
   vdt v_mhu_A1  = VEC_LDS(0, &mhu_a[7]);
   vdt v_mhu_A2  = VEC_LDS(0, &mhu_a[8]);
   vdt v_mhu_B1 = VEC_LDS(0, &mhu_a[14]);
   vdt v_mhu_B2 = VEC_LDS(0, &mhu_a[15]);
   vdt v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[11]);
   vdt v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[12]);
   vdt v_sum1a =  vec_madd(v_xa, v_mhu_A2, v_mhu_A1);
   vdt v_sum2a =  vec_madd(v_xa, v_mhu_B2, v_mhu_B1);
   vdt v_sum3a =  vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);
// 2
   v_mhu_A1  = VEC_LDS(0, &mhu_a[5]);
   v_mhu_A2  = VEC_LDS(0, &mhu_a[6]);
   v_mhu_B1 = VEC_LDS(0, &mhu_a[12]);
   v_mhu_B2 = VEC_LDS(0, &mhu_a[13]);
   v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[9]);
   v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[10]);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
// 3
   v_mhu_A1  = VEC_LDS(0, &mhu_a[3]);
   v_mhu_A2  = VEC_LDS(0, &mhu_a[4]);
   v_mhu_B1 = VEC_LDS(0, &mhu_a[10]);
   v_mhu_B2 = VEC_LDS(0, &mhu_a[11]);
   v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[7]);
   v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[8]);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
// 4
   v_mhu_A1  = VEC_LDS(0, &mhu_a[1]);
   v_mhu_A2  = VEC_LDS(0, &mhu_a[2]);
   v_mhu_B2  = VEC_LDS(0, &mhu_a[9]);
   v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[5]);
   v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[6]);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
// 5
   v_mhu_A2  = VEC_LDS(0, &mhu_a[0]);
   v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[3]);
   v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[4]);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
// 6
   v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[1]);
   v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[2]);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
// 7
   v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[0]);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);

// 1st
   vdt temp1a    = vec_swdiv_nochk(v_sum1a,v_sum2a); 
   vdt temp2a   = vec_ld(0, &g[ii]);  
   temp1a = vec_sub(temp1a, temp2a);  
   temp2a       = vec_madd(temp1a, v_sum3a, temp2a);
   vec_st(temp2a, 0, &g[ii]); //BODYEND
 }

 // the elements in the end of the vector should be ignored
 int leftover = nCells - ii;
 if (leftover) { 
      vdt v_xa   = vec_ld(0, &VM[ii]);  //BODYCOM1 
//BODYCOM1 
// 1//BODYCOM1 
   vdt v_mhu_A1  = VEC_LDS(0, &mhu_a[7]);//BODYCOM1 
   vdt v_mhu_A2  = VEC_LDS(0, &mhu_a[8]);//BODYCOM1 
   vdt v_mhu_B1 = VEC_LDS(0, &mhu_a[14]);//BODYCOM1 
   vdt v_mhu_B2 = VEC_LDS(0, &mhu_a[15]);//BODYCOM1 
   vdt v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[11]);//BODYCOM1 
   vdt v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[12]);//BODYCOM1 
   vdt v_sum1a =  vec_madd(v_xa, v_mhu_A2, v_mhu_A1);//BODYCOM1 
   vdt v_sum2a =  vec_madd(v_xa, v_mhu_B2, v_mhu_B1);//BODYCOM1 
   vdt v_sum3a =  vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);//BODYCOM1 
// 2//BODYCOM1 
   v_mhu_A1  = VEC_LDS(0, &mhu_a[5]);//BODYCOM1 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[6]);//BODYCOM1 
   v_mhu_B1 = VEC_LDS(0, &mhu_a[12]);//BODYCOM1 
   v_mhu_B2 = VEC_LDS(0, &mhu_a[13]);//BODYCOM1 
   v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[9]);//BODYCOM1 
   v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[10]);//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM1 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM1 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
// 3//BODYCOM1 
   v_mhu_A1  = VEC_LDS(0, &mhu_a[3]);//BODYCOM1 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[4]);//BODYCOM1 
   v_mhu_B1 = VEC_LDS(0, &mhu_a[10]);//BODYCOM1 
   v_mhu_B2 = VEC_LDS(0, &mhu_a[11]);//BODYCOM1 
   v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[7]);//BODYCOM1 
   v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[8]);//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM1 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM1 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
// 4//BODYCOM1 
   v_mhu_A1  = VEC_LDS(0, &mhu_a[1]);//BODYCOM1 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[2]);//BODYCOM1 
   v_mhu_B2  = VEC_LDS(0, &mhu_a[9]);//BODYCOM1 
   v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[5]);//BODYCOM1 
   v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[6]);//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM1 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
// 5//BODYCOM1 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[0]);//BODYCOM1 
   v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[3]);//BODYCOM1 
   v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[4]);//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
// 6//BODYCOM1 
   v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[1]);//BODYCOM1 
   v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[2]);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
// 7//BODYCOM1 
   v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[0]);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
//BODYCOM1 
// 1st//BODYCOM1 
   vdt temp1a    = vec_swdiv_nochk(v_sum1a,v_sum2a); //BODYCOM1 
   vdt temp2a   = vec_ld(0, &g[ii]);  //BODYCOM1 
   temp1a = vec_sub(temp1a, temp2a);  //BODYCOM1 
   temp2a       = vec_madd(temp1a, v_sum3a, temp2a);//BODYCOM1 

   
   //vec_st(temp2a, 0, &g[ii]);
   switch (leftover) {
   case 1:
     g[ii+0] = temp2a get [0];
     break;

   case 2:
     g[ii+0] = temp2a get [0];
     g[ii+1] = temp2a get [1];
     break;

   case 3:
     g[ii+0] = temp2a get [0];
     g[ii+1] = temp2a get [1];
     g[ii+2] = temp2a get [2];
     break;

     default:
    //println("BAD LEFTOVER: %d\n", leftover); fflush(stdout);
    assert(0);
   }
 }

}


void update_Xr1Gate_v1(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a)
{
  //update_Xr1Gate(dt, nCells, VM, g, mhu_a, tauR_a); return;
  typedef vector4double vdt;

int gateIndex=3;
int mhu_l= 5 ;
int mhu_m= 8;
int tauR_l= 1;
int tauR_m=13;

/*
 int mhu_k  = 12;  
 int tauR_k = 13;  
*/

 double  tauRdt_a[13];
         tauRdt_a[0]  = tauR_a[0]*dt,
         tauRdt_a[1]  = tauR_a[1]*dt,
         tauRdt_a[2]  = tauR_a[2]*dt,
         tauRdt_a[3]  = tauR_a[3]*dt,
         tauRdt_a[4]  = tauR_a[4]*dt,
         tauRdt_a[5]  = tauR_a[5]*dt,
         tauRdt_a[6]  = tauR_a[6]*dt, 
         tauRdt_a[7]  = tauR_a[7]*dt,
         tauRdt_a[8]  = tauR_a[8]*dt,
         tauRdt_a[9]  = tauR_a[9]*dt,
         tauRdt_a[10] = tauR_a[10]*dt,
         tauRdt_a[11] = tauR_a[11]*dt,
         tauRdt_a[12] = tauR_a[12]*dt;

  vdt v_mhu_a6 = VEC_LDS(0, &mhu_a[6]);
  vdt v_mhu_a7 = VEC_LDS(0, &mhu_a[7]);
  vdt v_sum4   = VEC_LDS(0, &tauR_a[13]);

  vdt v_TWO = vec_splats(2.0);

  int ii=0;

 // 0 aligned
 // 1: load produces (JUNK, &VM[0], &VM[1], &VM[2])
 // 2: load produces (JUNK, JUNK, &VM[0], &VM[1])
 // 3: load produces (JUNK, JUNK, JUNK, &VM[0])
 int misalignment = ((long int)&VM[0] % 32) / 8;
 assert(misalignment == 0);

 if (misalignment) { 
      vdt v_xa = vec_ld(0, &VM[ii]);  //BODYCOM0 
//BODYCOM0 
// 1//BODYCOM0 
 /*//BODYCOM0 
   sum1[ii] = mhu_a[7];//BODYCOM0 
   sum1[ii] = mhu_a[6]     + VM[ii]*sum1[ii];//BODYCOM0 
   sum2[ii] = mhu_a[12];//BODYCOM0 
   sum2[ii] = mhu_a[11]    + VM[ii]*sum2[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[12];//BODYCOM0 
   sum3[ii] = tauRdt_a[11] + VM[ii]*sum3[ii];//BODYCOM0 
*///BODYCOM0 
   vdt v_mhu_A1;//BODYCOM0 
   vdt v_mhu_A2;//BODYCOM0 
   vdt v_mhu_B1    = VEC_LDS(0, &mhu_a[11]);//BODYCOM0 
   vdt v_mhu_B2    = VEC_LDS(0, &mhu_a[12]);//BODYCOM0 
   vdt v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[11]);//BODYCOM0 
   vdt v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[12]);//BODYCOM0 
//BODYCOM0 
   vdt v_sum1a = vec_madd(v_xa, v_mhu_a7,    v_mhu_a6);//BODYCOM0 
   vdt v_sum2a = vec_madd(v_xa, v_mhu_B2,    v_mhu_B1);//BODYCOM0 
   vdt v_sum3a = vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);//BODYCOM0 
//BODYCOM0 
// 2  //BODYCOM0 
/*//BODYCOM0 
   sum1[ii] = mhu_a[5]     + VM[ii]*sum1[ii];//BODYCOM0 
   sum1[ii] = mhu_a[4]     + VM[ii]*sum1[ii];//BODYCOM0 
   sum2[ii] = mhu_a[10]    + VM[ii]*sum2[ii];//BODYCOM0 
   sum2[ii] = mhu_a[9]     + VM[ii]*sum2[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[10] + VM[ii]*sum3[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[9]  + VM[ii]*sum3[ii];//BODYCOM0 
*///BODYCOM0 
   v_mhu_A1    = VEC_LDS(0, &mhu_a[4]);//BODYCOM0 
   v_mhu_A2    = VEC_LDS(0, &mhu_a[5]);//BODYCOM0 
   v_mhu_B1    = VEC_LDS(0, &mhu_a[9]);//BODYCOM0 
   v_mhu_B2    = VEC_LDS(0, &mhu_a[10]);//BODYCOM0 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[9]);//BODYCOM0 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[10]);//BODYCOM0 
//BODYCOM0 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM0 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM0 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM0 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
//BODYCOM0 
// 3  //BODYCOM0 
/*//BODYCOM0 
   sum1[ii] = mhu_a[3]    + VM[ii]*sum1[ii];//BODYCOM0 
   sum1[ii] = mhu_a[2]    + VM[ii]*sum1[ii];//BODYCOM0 
   sum2[ii] = mhu_a[8]    + VM[ii]*sum2[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[8] + VM[ii]*sum3[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[7] + VM[ii]*sum3[ii];//BODYCOM0 
*///BODYCOM0 
   v_mhu_A1    = VEC_LDS(0, &mhu_a[2]);//BODYCOM0 
   v_mhu_A2    = VEC_LDS(0, &mhu_a[3]);//BODYCOM0 
   v_mhu_B1    = VEC_LDS(0, &mhu_a[8]);//BODYCOM0 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[7]);//BODYCOM0 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[8]);//BODYCOM0 
//BODYCOM0 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM0 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM0 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
//BODYCOM0 
// 4 //BODYCOM0 
/*//BODYCOM0 
   sum1[ii] = mhu_a[1]    + VM[ii]*sum1[ii];//BODYCOM0 
   sum1[ii] = mhu_a[0]    + VM[ii]*sum1[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[6] + VM[ii]*sum3[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[5] + VM[ii]*sum3[ii];//BODYCOM0 
*///BODYCOM0 
   v_mhu_A1    = VEC_LDS(0, &mhu_a[0]);//BODYCOM0 
   v_mhu_A2    = VEC_LDS(0, &mhu_a[1]);//BODYCOM0 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[5]);//BODYCOM0 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[6]);//BODYCOM0 
//BODYCOM0 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM0 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
//BODYCOM0 
// 5 //BODYCOM0 
/*//BODYCOM0 
   sum3[ii] = tauRdt_a[4] + VM[ii]*sum3[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[3] + VM[ii]*sum3[ii];//BODYCOM0 
*///BODYCOM0 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[3]);//BODYCOM0 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[4]);//BODYCOM0 
//BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
//BODYCOM0 
// 6 //BODYCOM0 
/*//BODYCOM0 
   sum3[ii] = tauRdt_a[2] + VM[ii]*sum3[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[1] + VM[ii]*sum3[ii];//BODYCOM0 
*///BODYCOM0 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[1]);//BODYCOM0 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[2]);//BODYCOM0 
//BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
//BODYCOM0 
// 7 //BODYCOM0 
/*//BODYCOM0 
   sum3[ii]  =  tauRdt_a[0]  + VM[ii]*sum3[ii];//BODYCOM0 
*///BODYCOM0 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[0]);//BODYCOM0 
//BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
//BODYCOM0 
/*//BODYCOM0 
#if (1) //BODYCOM0 
    double tauRdt= sum3/sum4; //BODYCOM0 
    double mhu= sum1/sum2;//BODYCOM0 
    g[ii] +=  mhu*tauRdt - g[ii]*tauRdt;//BODYCOM0 
#else//BODYCOM0 
   double x = sum2*sum4; //BODYCOM0 
   double f = recipApprox(x); //BODYCOM0 
   double a = sum1-g[ii]*sum2; //BODYCOM0 
   double b = sum3*f*(2.0-f *x) ;//BODYCOM0 
   g[ii] +=  a*b; //BODYCOM0 
#endif//BODYCOM0 
*///BODYCOM0 
//BODYCOM0 
//BODYCOM0 
// Still need to reuse the register I have above .//BODYCOM0 
// A//BODYCOM0 
   v_xa          = vec_swdiv_nochk(v_sum3a, v_sum4); //BODYCOM0 
   v_mhu_A2      = vec_swdiv_nochk(v_sum1a, v_sum2a); //BODYCOM0 
   v_mhu_B2      = vec_ld(0, &g[ii]);  //BODYCOM0 
   vdt v_tauR_D2 = vec_sub(v_mhu_A2, v_mhu_B2);  //BODYCOM0 
   v_mhu_B2      = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);//BODYCOM0 

   
   //vec_st(v_mhu_B2, 0, &g[ii]);
   switch (misalignment) {
   case 1:
     g[ii+0] = v_mhu_B2 get [1];
     g[ii+1] = v_mhu_B2 get [2];
     g[ii+2] = v_mhu_B2 get [3];
     break;

   case 2:
     g[ii+0] = v_mhu_B2 get [2];
     g[ii+1] = v_mhu_B2 get [3];
     break;

   case 3:
     g[ii+0] = v_mhu_B2 get [3];
     break;

     default:
    //println("BAD MISALIGNMENT: %d\n", misalignment); fflush(stdout);
     assert(0);
   }

   ii = 4 - misalignment;
 }

 for (ii;ii<TTROUND(nCells,8);ii+=8)
 { 
   vdt v_xa = vec_ld(0, &VM[ii]);  
   vdt v_xb = vec_ld(0, &VM[ii+4]);  
// 1
 /*
   sum1[ii] = mhu_a[7];
   sum1[ii] = mhu_a[6]     + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[12];
   sum2[ii] = mhu_a[11]    + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[12];
   sum3[ii] = tauRdt_a[11] + VM[ii]*sum3[ii];
*/
   vdt v_mhu_A1;
   vdt v_mhu_A2;
   vdt v_mhu_B1    = VEC_LDS(0, &mhu_a[11]);
   vdt v_mhu_B2    = VEC_LDS(0, &mhu_a[12]);
   vdt v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[11]);
   vdt v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[12]);

   vdt v_sum1a = vec_madd(v_xa, v_mhu_a7,    v_mhu_a6);
   vdt v_sum1b = vec_madd(v_xb, v_mhu_a7,    v_mhu_a6);
   vdt v_sum2a = vec_madd(v_xa, v_mhu_B2,    v_mhu_B1);
   vdt v_sum2b = vec_madd(v_xb, v_mhu_B2,    v_mhu_B1);
   vdt v_sum3a = vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);
   vdt v_sum3b = vec_madd(v_xb, v_tauRdt_C2, v_tauRdt_C1);

// 2  
/*
   sum1[ii] = mhu_a[5]     + VM[ii]*sum1[ii];
   sum1[ii] = mhu_a[4]     + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[10]    + VM[ii]*sum2[ii];
   sum2[ii] = mhu_a[9]     + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[10] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[9]  + VM[ii]*sum3[ii];
*/
   v_mhu_A1    = VEC_LDS(0, &mhu_a[4]);
   v_mhu_A2    = VEC_LDS(0, &mhu_a[5]);
   v_mhu_B1    = VEC_LDS(0, &mhu_a[9]);
   v_mhu_B2    = VEC_LDS(0, &mhu_a[10]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[9]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[10]);

   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);

   v_sum1b = vec_madd(v_xb, v_sum1b, v_mhu_A2);
   v_sum1b = vec_madd(v_xb, v_sum1b, v_mhu_A1);
   v_sum2b = vec_madd(v_xb, v_sum2b, v_mhu_B2);
   v_sum2b = vec_madd(v_xb, v_sum2b, v_mhu_B1); 
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C1);

// 3  
/*
   sum1[ii] = mhu_a[3]    + VM[ii]*sum1[ii];
   sum1[ii] = mhu_a[2]    + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[8]    + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[8] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[7] + VM[ii]*sum3[ii];
*/
   v_mhu_A1    = VEC_LDS(0, &mhu_a[2]);
   v_mhu_A2    = VEC_LDS(0, &mhu_a[3]);
   v_mhu_B1    = VEC_LDS(0, &mhu_a[8]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[7]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[8]);

   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);

   v_sum1b = vec_madd(v_xb, v_sum1b, v_mhu_A2);
   v_sum1b = vec_madd(v_xb, v_sum1b, v_mhu_A1);
   v_sum2b = vec_madd(v_xb, v_sum2b, v_mhu_B1); 
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C1);

// 4 
/*
   sum1[ii] = mhu_a[1]    + VM[ii]*sum1[ii];
   sum1[ii] = mhu_a[0]    + VM[ii]*sum1[ii];
   sum3[ii] = tauRdt_a[6] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[5] + VM[ii]*sum3[ii];
*/
   v_mhu_A1    = VEC_LDS(0, &mhu_a[0]);
   v_mhu_A2    = VEC_LDS(0, &mhu_a[1]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[5]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[6]);

   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);

   v_sum1b = vec_madd(v_xb, v_sum1b, v_mhu_A2);
   v_sum1b = vec_madd(v_xb, v_sum1b, v_mhu_A1);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C1);

// 5 
/*
   sum3[ii] = tauRdt_a[4] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[3] + VM[ii]*sum3[ii];
*/
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[3]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[4]);

   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);

   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C1);

// 6 
/*
   sum3[ii] = tauRdt_a[2] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[1] + VM[ii]*sum3[ii];
*/
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[1]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[2]);

   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);

   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C1);

// 7 
/*
   sum3[ii]  =  tauRdt_a[0]  + VM[ii]*sum3[ii];
*/
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[0]);

   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);

   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C1);

/*
#if (1) 
    double tauRdt= sum3/sum4; 
    double mhu= sum1/sum2;
    g[ii] +=  mhu*tauRdt - g[ii]*tauRdt;
#else
   double x = sum2*sum4; 
   double f = recipApprox(x); 
   double a = sum1-g[ii]*sum2; 
   double b = sum3*f*(2.0-f *x) ;
   g[ii] +=  a*b; 
#endif
*/

#if (1) 
// Still need to reuse the register I have above .
// A
   v_xa          = vec_swdiv_nochk(v_sum3a, v_sum4); 
   v_mhu_A2      = vec_swdiv_nochk(v_sum1a, v_sum2a); 
   v_mhu_B2      = vec_ld(0, &g[ii]);  
   vdt v_tauR_D2 = vec_sub(v_mhu_A2, v_mhu_B2);  
   v_mhu_B2      = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);
   vec_st(v_mhu_B2, 0, &g[ii]);

// B
   v_xb          = vec_swdiv_nochk(v_sum3b, v_sum4); 
   v_mhu_A1      = vec_swdiv_nochk(v_sum1b, v_sum2b); 
   v_mhu_B1      = vec_ld(0, &g[ii+4]);  
   vdt v_tauR_D1 = vec_sub(v_mhu_A1, v_mhu_B1);  
   v_mhu_B1  = vec_madd(v_tauR_D1, v_xb, v_mhu_B1);
   vec_st(v_mhu_B1, 0, &g[ii+4]);
#else
//***************** NEED TO DEBUG
// A
   v_tauRdt_C2 = vec_ld(0, &g[ii]); 
 
   v_xa      = vec_mul(v_sum2a, v_sum4);
   v_mhu_A2  = vec_re(v_xa); 
   v_mhu_B2  = vec_nmsub(v_tauRdt_C2, v_sum2a, v_sum1a);
   vdt v_tauR_D2 = vec_nmsub(v_mhu_A2, v_xa, v_TWO);
   v_tauR_D2 = vec_mul(v_mhu_A2,v_tauR_D2 );
   v_tauR_D2 = vec_mul(v_sum3a, v_tauR_D2);

   v_tauRdt_C2  = vec_madd(v_mhu_B2,v_tauR_D2 , v_tauRdt_C2);

   vec_st(v_tauRdt_C2, 0, &g[ii]);
// B 
   v_tauRdt_C1 = vec_ld(0, &g[ii+4]); 

   v_xb      = vec_mul(v_sum2b, v_sum4);
   v_mhu_A1  = vec_re(v_xb); 
   v_mhu_B1  = vec_nmsub(v_tauRdt_C1, v_sum2b, v_sum1b);
   vdt v_tauR_D1 = vec_nmsub(v_mhu_A1, v_xb, v_TWO);
   v_tauR_D1 = vec_mul(v_mhu_A1,v_tauR_D1 );
   v_tauR_D1 = vec_mul(v_sum3b, v_tauR_D1);

   v_tauRdt_C1  = vec_madd(v_mhu_B1,v_tauR_D1 , v_tauRdt_C1);
   vec_st(v_tauRdt_C1, 0, &g[ii+4]);

#endif

  } 
 for (ii;ii<TTROUND(nCells,4);ii+=4)
 { 
  //BODYSTART
   vdt v_xa = vec_ld(0, &VM[ii]);  

// 1
 /*
   sum1[ii] = mhu_a[7];
   sum1[ii] = mhu_a[6]     + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[12];
   sum2[ii] = mhu_a[11]    + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[12];
   sum3[ii] = tauRdt_a[11] + VM[ii]*sum3[ii];
*/
   vdt v_mhu_A1;
   vdt v_mhu_A2;
   vdt v_mhu_B1    = VEC_LDS(0, &mhu_a[11]);
   vdt v_mhu_B2    = VEC_LDS(0, &mhu_a[12]);
   vdt v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[11]);
   vdt v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[12]);

   vdt v_sum1a = vec_madd(v_xa, v_mhu_a7,    v_mhu_a6);
   vdt v_sum2a = vec_madd(v_xa, v_mhu_B2,    v_mhu_B1);
   vdt v_sum3a = vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);

// 2  
/*
   sum1[ii] = mhu_a[5]     + VM[ii]*sum1[ii];
   sum1[ii] = mhu_a[4]     + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[10]    + VM[ii]*sum2[ii];
   sum2[ii] = mhu_a[9]     + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[10] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[9]  + VM[ii]*sum3[ii];
*/
   v_mhu_A1    = VEC_LDS(0, &mhu_a[4]);
   v_mhu_A2    = VEC_LDS(0, &mhu_a[5]);
   v_mhu_B1    = VEC_LDS(0, &mhu_a[9]);
   v_mhu_B2    = VEC_LDS(0, &mhu_a[10]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[9]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[10]);

   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);

// 3  
/*
   sum1[ii] = mhu_a[3]    + VM[ii]*sum1[ii];
   sum1[ii] = mhu_a[2]    + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[8]    + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[8] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[7] + VM[ii]*sum3[ii];
*/
   v_mhu_A1    = VEC_LDS(0, &mhu_a[2]);
   v_mhu_A2    = VEC_LDS(0, &mhu_a[3]);
   v_mhu_B1    = VEC_LDS(0, &mhu_a[8]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[7]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[8]);

   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);

// 4 
/*
   sum1[ii] = mhu_a[1]    + VM[ii]*sum1[ii];
   sum1[ii] = mhu_a[0]    + VM[ii]*sum1[ii];
   sum3[ii] = tauRdt_a[6] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[5] + VM[ii]*sum3[ii];
*/
   v_mhu_A1    = VEC_LDS(0, &mhu_a[0]);
   v_mhu_A2    = VEC_LDS(0, &mhu_a[1]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[5]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[6]);

   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);

// 5 
/*
   sum3[ii] = tauRdt_a[4] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[3] + VM[ii]*sum3[ii];
*/
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[3]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[4]);

   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);

// 6 
/*
   sum3[ii] = tauRdt_a[2] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[1] + VM[ii]*sum3[ii];
*/
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[1]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[2]);

   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);

// 7 
/*
   sum3[ii]  =  tauRdt_a[0]  + VM[ii]*sum3[ii];
*/
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[0]);

   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);

/*
#if (1) 
    double tauRdt= sum3/sum4; 
    double mhu= sum1/sum2;
    g[ii] +=  mhu*tauRdt - g[ii]*tauRdt;
#else
   double x = sum2*sum4; 
   double f = recipApprox(x); 
   double a = sum1-g[ii]*sum2; 
   double b = sum3*f*(2.0-f *x) ;
   g[ii] +=  a*b; 
#endif
*/


// Still need to reuse the register I have above .
// A
   v_xa          = vec_swdiv_nochk(v_sum3a, v_sum4); 
   v_mhu_A2      = vec_swdiv_nochk(v_sum1a, v_sum2a); 
   v_mhu_B2      = vec_ld(0, &g[ii]);  
   vdt v_tauR_D2 = vec_sub(v_mhu_A2, v_mhu_B2);  
   v_mhu_B2      = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);
   vec_st(v_mhu_B2, 0, &g[ii]); //BODYEND

#if 0
//***************** NEED TO DEBUG
// A
   v_tauRdt_C2 = vec_ld(0, &g[ii]); 
 
   v_xa      = vec_mul(v_sum2a, v_sum4);
   v_mhu_A2  = vec_re(v_xa); 
   v_mhu_B2  = vec_nmsub(v_tauRdt_C2, v_sum2a, v_sum1a);
   vdt v_tauR_D2 = vec_nmsub(v_mhu_A2, v_xa, v_TWO);
   v_tauR_D2 = vec_mul(v_mhu_A2,v_tauR_D2 );
   v_tauR_D2 = vec_mul(v_sum3a, v_tauR_D2);

   v_tauRdt_C2  = vec_madd(v_mhu_B2,v_tauR_D2 , v_tauRdt_C2);

   vec_st(v_tauRdt_C2, 0, &g[ii]);
#endif

  }

 // the elements in the end of the vector should be ignored
 int leftover = nCells - ii;
 if (leftover) { 
      vdt v_xa = vec_ld(0, &VM[ii]);  //BODYCOM1 
//BODYCOM1 
// 1//BODYCOM1 
 /*//BODYCOM1 
   sum1[ii] = mhu_a[7];//BODYCOM1 
   sum1[ii] = mhu_a[6]     + VM[ii]*sum1[ii];//BODYCOM1 
   sum2[ii] = mhu_a[12];//BODYCOM1 
   sum2[ii] = mhu_a[11]    + VM[ii]*sum2[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[12];//BODYCOM1 
   sum3[ii] = tauRdt_a[11] + VM[ii]*sum3[ii];//BODYCOM1 
*///BODYCOM1 
   vdt v_mhu_A1;//BODYCOM1 
   vdt v_mhu_A2;//BODYCOM1 
   vdt v_mhu_B1    = VEC_LDS(0, &mhu_a[11]);//BODYCOM1 
   vdt v_mhu_B2    = VEC_LDS(0, &mhu_a[12]);//BODYCOM1 
   vdt v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[11]);//BODYCOM1 
   vdt v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[12]);//BODYCOM1 
//BODYCOM1 
   vdt v_sum1a = vec_madd(v_xa, v_mhu_a7,    v_mhu_a6);//BODYCOM1 
   vdt v_sum2a = vec_madd(v_xa, v_mhu_B2,    v_mhu_B1);//BODYCOM1 
   vdt v_sum3a = vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);//BODYCOM1 
//BODYCOM1 
// 2  //BODYCOM1 
/*//BODYCOM1 
   sum1[ii] = mhu_a[5]     + VM[ii]*sum1[ii];//BODYCOM1 
   sum1[ii] = mhu_a[4]     + VM[ii]*sum1[ii];//BODYCOM1 
   sum2[ii] = mhu_a[10]    + VM[ii]*sum2[ii];//BODYCOM1 
   sum2[ii] = mhu_a[9]     + VM[ii]*sum2[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[10] + VM[ii]*sum3[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[9]  + VM[ii]*sum3[ii];//BODYCOM1 
*///BODYCOM1 
   v_mhu_A1    = VEC_LDS(0, &mhu_a[4]);//BODYCOM1 
   v_mhu_A2    = VEC_LDS(0, &mhu_a[5]);//BODYCOM1 
   v_mhu_B1    = VEC_LDS(0, &mhu_a[9]);//BODYCOM1 
   v_mhu_B2    = VEC_LDS(0, &mhu_a[10]);//BODYCOM1 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[9]);//BODYCOM1 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[10]);//BODYCOM1 
//BODYCOM1 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM1 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM1 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM1 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
//BODYCOM1 
// 3  //BODYCOM1 
/*//BODYCOM1 
   sum1[ii] = mhu_a[3]    + VM[ii]*sum1[ii];//BODYCOM1 
   sum1[ii] = mhu_a[2]    + VM[ii]*sum1[ii];//BODYCOM1 
   sum2[ii] = mhu_a[8]    + VM[ii]*sum2[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[8] + VM[ii]*sum3[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[7] + VM[ii]*sum3[ii];//BODYCOM1 
*///BODYCOM1 
   v_mhu_A1    = VEC_LDS(0, &mhu_a[2]);//BODYCOM1 
   v_mhu_A2    = VEC_LDS(0, &mhu_a[3]);//BODYCOM1 
   v_mhu_B1    = VEC_LDS(0, &mhu_a[8]);//BODYCOM1 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[7]);//BODYCOM1 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[8]);//BODYCOM1 
//BODYCOM1 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM1 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM1 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
//BODYCOM1 
// 4 //BODYCOM1 
/*//BODYCOM1 
   sum1[ii] = mhu_a[1]    + VM[ii]*sum1[ii];//BODYCOM1 
   sum1[ii] = mhu_a[0]    + VM[ii]*sum1[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[6] + VM[ii]*sum3[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[5] + VM[ii]*sum3[ii];//BODYCOM1 
*///BODYCOM1 
   v_mhu_A1    = VEC_LDS(0, &mhu_a[0]);//BODYCOM1 
   v_mhu_A2    = VEC_LDS(0, &mhu_a[1]);//BODYCOM1 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[5]);//BODYCOM1 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[6]);//BODYCOM1 
//BODYCOM1 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM1 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
//BODYCOM1 
// 5 //BODYCOM1 
/*//BODYCOM1 
   sum3[ii] = tauRdt_a[4] + VM[ii]*sum3[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[3] + VM[ii]*sum3[ii];//BODYCOM1 
*///BODYCOM1 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[3]);//BODYCOM1 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[4]);//BODYCOM1 
//BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
//BODYCOM1 
// 6 //BODYCOM1 
/*//BODYCOM1 
   sum3[ii] = tauRdt_a[2] + VM[ii]*sum3[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[1] + VM[ii]*sum3[ii];//BODYCOM1 
*///BODYCOM1 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[1]);//BODYCOM1 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[2]);//BODYCOM1 
//BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
//BODYCOM1 
// 7 //BODYCOM1 
/*//BODYCOM1 
   sum3[ii]  =  tauRdt_a[0]  + VM[ii]*sum3[ii];//BODYCOM1 
*///BODYCOM1 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[0]);//BODYCOM1 
//BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
//BODYCOM1 
/*//BODYCOM1 
#if (1) //BODYCOM1 
    double tauRdt= sum3/sum4; //BODYCOM1 
    double mhu= sum1/sum2;//BODYCOM1 
    g[ii] +=  mhu*tauRdt - g[ii]*tauRdt;//BODYCOM1 
#else//BODYCOM1 
   double x = sum2*sum4; //BODYCOM1 
   double f = recipApprox(x); //BODYCOM1 
   double a = sum1-g[ii]*sum2; //BODYCOM1 
   double b = sum3*f*(2.0-f *x) ;//BODYCOM1 
   g[ii] +=  a*b; //BODYCOM1 
#endif//BODYCOM1 
*///BODYCOM1 
//BODYCOM1 
//BODYCOM1 
// Still need to reuse the register I have above .//BODYCOM1 
// A//BODYCOM1 
   v_xa          = vec_swdiv_nochk(v_sum3a, v_sum4); //BODYCOM1 
   v_mhu_A2      = vec_swdiv_nochk(v_sum1a, v_sum2a); //BODYCOM1 
   v_mhu_B2      = vec_ld(0, &g[ii]);  //BODYCOM1 
   vdt v_tauR_D2 = vec_sub(v_mhu_A2, v_mhu_B2);  //BODYCOM1 
   v_mhu_B2      = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);//BODYCOM1 

   
   //vec_st(v_mhu_B2, 0, &g[ii]);
   //printf("%d: leftover=%d\n",gateIndex,leftover); fflush(stdout); 
   switch (leftover) {
   case 1:
     g[ii+0] = v_mhu_B2 get [0];
     break;

   case 2:
     g[ii+0] = v_mhu_B2 get [0];
     g[ii+1] = v_mhu_B2 get [1];
     break;

   case 3:
     g[ii+0] = v_mhu_B2 get [0];
     g[ii+1] = v_mhu_B2 get [1];
     g[ii+2] = v_mhu_B2 get [2];
     break;

     default:
    //println("BAD LEFTOVER: %d\n", leftover); fflush(stdout);
    assert(0);
   }
 }

}


void update_Xr2Gate_v1(double dt, int nCells,  double *VM, double *g, double *mhu_a, double *tauR_a)
{
  //update_Xr2Gate(dt, nCells,  VM, g, mhu_a, tauR_a); return;
  typedef vector4double vdt;

int gateIndex=4;
int mhu_l=1 ;
int mhu_m= 10;
int tauR_l=1;
int tauR_m=10;
   
 double  tauRdt_a[10];



  vdt v_mhu_a0  = VEC_LDS(0, &mhu_a[0]);
  vdt v_mhu_a1  = VEC_LDS(0, &mhu_a[1]);
  vdt v_mhu_a2  = VEC_LDS(0, &mhu_a[2]);
  vdt v_mhu_a3  = VEC_LDS(0, &mhu_a[3]);
  vdt v_mhu_a4  = VEC_LDS(0, &mhu_a[4]);
  vdt v_mhu_a5  = VEC_LDS(0, &mhu_a[5]);
  vdt v_mhu_a6  = VEC_LDS(0, &mhu_a[6]);
  vdt v_mhu_a7  = VEC_LDS(0, &mhu_a[7]);
  vdt v_mhu_a8  = VEC_LDS(0, &mhu_a[8]);
  vdt v_mhu_a9  = VEC_LDS(0, &mhu_a[9]);


double   tauRdt_a0  = tauR_a[0]*dt,
         tauRdt_a1  = tauR_a[1]*dt,
         tauRdt_a2  = tauR_a[2]*dt,
         tauRdt_a3  = tauR_a[3]*dt,
         tauRdt_a4  = tauR_a[4]*dt,
         tauRdt_a5  = tauR_a[5]*dt,
         tauRdt_a6  = tauR_a[6]*dt, 
         tauRdt_a7  = tauR_a[7]*dt,
         tauRdt_a8  = tauR_a[8]*dt,
         tauRdt_a9  = tauR_a[9]*dt;

  vdt v_tauRdt_a0   = vec_splats(tauRdt_a0);
  vdt v_tauRdt_a1   = vec_splats(tauRdt_a1);
  vdt v_tauRdt_a2   = vec_splats(tauRdt_a2);
  vdt v_tauRdt_a3   = vec_splats(tauRdt_a3);
  vdt v_tauRdt_a4   = vec_splats(tauRdt_a4);
  vdt v_tauRdt_a5   = vec_splats(tauRdt_a5);
  vdt v_tauRdt_a6   = vec_splats(tauRdt_a6);
  vdt v_tauRdt_a7   = vec_splats(tauRdt_a7);
  vdt v_tauRdt_a8   = vec_splats(tauRdt_a8);
  vdt v_tauRdt_a9   = vec_splats(tauRdt_a9);


  // vdt v_TWO = vec_splats(2.0);

 int mhu_k  = 10;  
 int tauR_k = 10;  
 int ii=0;

 // 0 aligned
 // 1: load produces (JUNK, &VM[0], &VM[1], &VM[2])
 // 2: load produces (JUNK, JUNK, &VM[0], &VM[1])
 // 3: load produces (JUNK, JUNK, JUNK, &VM[0])
 int misalignment = ((long int)&VM[0] % 32) / 8;
 assert(misalignment == 0);

 if (misalignment) { 
      vdt v_v = vec_ld(0, &VM[ii]); //BODYCOM0 
//BODYCOM0 
   vdt v_sum1   =  vec_madd(v_v, v_mhu_a9, v_mhu_a8);//BODYCOM0 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a7);//BODYCOM0 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a6);//BODYCOM0 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a5);//BODYCOM0 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a4); //BODYCOM0 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a3);//BODYCOM0 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a2);//BODYCOM0 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a1);//BODYCOM0 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a0);//BODYCOM0 
//BODYCOM0 
   vdt v_sum3   =  vec_madd(v_v, v_tauRdt_a9, v_tauRdt_a8);//BODYCOM0 
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a7);//BODYCOM0 
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a6);//BODYCOM0 
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a5);//BODYCOM0 
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a4);//BODYCOM0 
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a3);//BODYCOM0 
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a2);//BODYCOM0 
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a1);//BODYCOM0 
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a0);//BODYCOM0 
//BODYCOM0 
   vdt v_g     = vec_ld(0, &g[ii]);  //BODYCOM0 
   vdt temp1     = vec_madd(v_sum1,   v_sum3,   v_g);//BODYCOM0 
   v_g  = vec_nmsub( v_g,   v_sum3,   temp1);//BODYCOM0 

   
   //vec_st(v_g, 0, &g[ii]);
   switch (misalignment) {
   case 1:
     g[ii+0] = v_g get [1];
     g[ii+1] = v_g get [2];
     g[ii+2] = v_g get [3];
     break;

   case 2:
     g[ii+0] = v_g get [2];
     g[ii+1] = v_g get [3];
     break;

   case 3:
     g[ii+0] = v_g get [3];
     break;

     default:
    //println("BAD MISALIGNMENT: %d\n", misalignment); fflush(stdout);
     assert(0);
   }

   ii = 4 - misalignment;
 }

 for (ii;ii<TTROUND(nCells,8);ii+=8)
 { 
   vdt v_v = vec_ld(0, &VM[ii]); 
   vdt v_v_4 = vec_ld(0, &VM[ii+4]);    

   vdt v_sum1   =  vec_madd(v_v, v_mhu_a9, v_mhu_a8);
   vdt v_sum1_4 =  vec_madd(v_v_4, v_mhu_a9, v_mhu_a8);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a7);
   v_sum1_4 =  vec_madd(v_v_4, v_sum1_4, v_mhu_a7);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a6);
   v_sum1_4 =  vec_madd(v_v_4, v_sum1_4, v_mhu_a6);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a5);
   v_sum1_4 =  vec_madd(v_v_4, v_sum1_4, v_mhu_a5);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a4); 
   v_sum1_4 =  vec_madd(v_v_4, v_sum1_4, v_mhu_a4); 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a3);
   v_sum1_4 =  vec_madd(v_v_4, v_sum1_4, v_mhu_a3);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a2);
   v_sum1_4 =  vec_madd(v_v_4, v_sum1_4, v_mhu_a2);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a1);
   v_sum1_4 =  vec_madd(v_v_4, v_sum1_4, v_mhu_a1);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a0);
   v_sum1_4 =  vec_madd(v_v_4, v_sum1_4, v_mhu_a0);

   vdt v_sum3   =  vec_madd(v_v, v_tauRdt_a9, v_tauRdt_a8);
   vdt v_sum3_4 =  vec_madd(v_v_4, v_tauRdt_a9, v_tauRdt_a8);
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a7);
   v_sum3_4 =  vec_madd(v_v_4, v_sum3_4, v_tauRdt_a7);
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a6);
   v_sum3_4 =  vec_madd(v_v_4, v_sum3_4, v_tauRdt_a6);
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a5);
   v_sum3_4 =  vec_madd(v_v_4, v_sum3_4, v_tauRdt_a5);
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a4);
   v_sum3_4 =  vec_madd(v_v_4, v_sum3_4, v_tauRdt_a4);
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a3);
   v_sum3_4 =  vec_madd(v_v_4, v_sum3_4, v_tauRdt_a3);
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a2);
   v_sum3_4 =  vec_madd(v_v_4, v_sum3_4, v_tauRdt_a2);
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a1);
   v_sum3_4 =  vec_madd(v_v_4, v_sum3_4, v_tauRdt_a1);
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a0);
   v_sum3_4 =  vec_madd(v_v_4, v_sum3_4, v_tauRdt_a0);

   vdt v_g     = vec_ld(0, &g[ii]);  
   vdt v_g_4   = vec_ld(0, &g[ii+4]);
   vdt temp1     = vec_madd(v_sum1,   v_sum3,   v_g);
   vdt temp1_4   = vec_madd(v_sum1_4, v_sum3_4, v_g_4);
   v_g   = vec_nmsub( v_g,   v_sum3,   temp1);
   v_g_4 = vec_nmsub( v_g_4, v_sum3_4, temp1_4);  
   vec_st(v_g, 0, &g[ii]);
   vec_st(v_g_4, 0, &g[ii+4]);
 }
 for (ii;ii<TTROUND(nCells,4);ii+=4)
 { 
   //BODYSTART
   vdt v_v = vec_ld(0, &VM[ii]); 

   vdt v_sum1   =  vec_madd(v_v, v_mhu_a9, v_mhu_a8);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a7);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a6);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a5);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a4); 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a3);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a2);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a1);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a0);

   vdt v_sum3   =  vec_madd(v_v, v_tauRdt_a9, v_tauRdt_a8);
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a7);
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a6);
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a5);
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a4);
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a3);
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a2);
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a1);
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a0);

   vdt v_g     = vec_ld(0, &g[ii]);  
   vdt temp1     = vec_madd(v_sum1,   v_sum3,   v_g);
   v_g  = vec_nmsub( v_g,   v_sum3,   temp1);
   vec_st(v_g, 0, &g[ii]); //BODYEND
 }

 // the elements in the end of the vector should be ignored
 int leftover = nCells - ii;
 if (leftover) { 
      vdt v_v = vec_ld(0, &VM[ii]); //BODYCOM1 
//BODYCOM1 
   vdt v_sum1   =  vec_madd(v_v, v_mhu_a9, v_mhu_a8);//BODYCOM1 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a7);//BODYCOM1 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a6);//BODYCOM1 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a5);//BODYCOM1 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a4); //BODYCOM1 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a3);//BODYCOM1 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a2);//BODYCOM1 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a1);//BODYCOM1 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a0);//BODYCOM1 
//BODYCOM1 
   vdt v_sum3   =  vec_madd(v_v, v_tauRdt_a9, v_tauRdt_a8);//BODYCOM1 
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a7);//BODYCOM1 
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a6);//BODYCOM1 
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a5);//BODYCOM1 
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a4);//BODYCOM1 
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a3);//BODYCOM1 
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a2);//BODYCOM1 
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a1);//BODYCOM1 
   v_sum3   =  vec_madd(v_v, v_sum3, v_tauRdt_a0);//BODYCOM1 
//BODYCOM1 
   vdt v_g     = vec_ld(0, &g[ii]);  //BODYCOM1 
   vdt temp1     = vec_madd(v_sum1,   v_sum3,   v_g);//BODYCOM1 
   v_g  = vec_nmsub( v_g,   v_sum3,   temp1);//BODYCOM1 

   
   //vec_st(v_g, 0, &g[ii]);
   switch (leftover) {
   case 1:
     g[ii+0] = v_g get [0];
     break;

   case 2:
     g[ii+0] = v_g get [0];
     g[ii+1] = v_g get [1];
     break;

   case 3:
     g[ii+0] = v_g get [0];
     g[ii+1] = v_g get [1];
     g[ii+2] = v_g get [2];
     break;

     default:
    //println("BAD LEFTOVER: %d\n", leftover); fflush(stdout);
    assert(0);
   }
 }

}


void update_XsGate_v1(double dt, int nCells,  double *VM, double *g, double *mhu_a, double *tauR_a)
{
  //update_XsGate(dt, nCells,  VM, g, mhu_a, tauR_a); return;

  typedef vector4double vdt;

int gateIndex= 5;
int mhu_l= 5 ;
int mhu_m= 5;
int tauR_l= 6;
int tauR_m= 9;

/*
 int mhu_k  = 9;  
 int tauR_k = 14  
*/
   
 double  tauRdt_a[9];
 // for (int j=8;j>=0;j--)    tauRdt_a[j] = tauR_a[j]*dt;
         tauRdt_a[0]  = tauR_a[0]*dt,
         tauRdt_a[1]  = tauR_a[1]*dt,
         tauRdt_a[2]  = tauR_a[2]*dt,
         tauRdt_a[3]  = tauR_a[3]*dt,
         tauRdt_a[4]  = tauR_a[4]*dt,
         tauRdt_a[5]  = tauR_a[5]*dt,
         tauRdt_a[6]  = tauR_a[6]*dt, 
         tauRdt_a[7]  = tauR_a[7]*dt,
         tauRdt_a[8]  = tauR_a[8]*dt;


 int ii=0;

 // 0 aligned
 // 1: load produces (JUNK, &VM[0], &VM[1], &VM[2])
 // 2: load produces (JUNK, JUNK, &VM[0], &VM[1])
 // 3: load produces (JUNK, JUNK, JUNK, &VM[0])
 int misalignment = ((long int)&VM[0] % 32) / 8;
 assert(misalignment == 0);

 if (misalignment) { 
      vdt v_xa   = vec_ld(0, &VM[ii]);  //BODYCOM0 
 //BODYCOM0 
// 1//BODYCOM0 
 /*//BODYCOM0 
   sum1[ii] =  mhu_a[4];//BODYCOM0 
   sum1[ii] =  mhu_a[3]    + VM[ii]*sum1[ii];//BODYCOM0 
   sum2[ii] =  mhu_a[9];//BODYCOM0 
   sum2[ii] =  mhu_a[8]    + VM[ii]*sum2[ii];//BODYCOM0 
   sum3[ii] =  tauRdt_a[8];//BODYCOM0 
   sum3[ii] =  tauRdt_a[7] + VM[ii]*sum3[ii];//BODYCOM0 
   sum4[ii] =  tauR_a[14]; //BODYCOM0 
   sum4[ii] =  tauR_a[13]  + VM[ii]*sum4[ii]; //BODYCOM0 
*///BODYCOM0 
  vdt v_mhu_A1    = VEC_LDS(0, &mhu_a[3]);//BODYCOM0 
  vdt v_mhu_A2    = VEC_LDS(0, &mhu_a[4]);//BODYCOM0 
  vdt v_mhu_B1    = VEC_LDS(0, &mhu_a[8]);//BODYCOM0 
  vdt v_mhu_B2    = VEC_LDS(0, &mhu_a[9]);//BODYCOM0 
  vdt v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[7]);//BODYCOM0 
  vdt v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[8]);//BODYCOM0 
  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[13]);//BODYCOM0 
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[14]);//BODYCOM0 
//BODYCOM0 
  vdt v_sum1a = vec_madd(v_xa, v_mhu_A2, v_mhu_A1);//BODYCOM0 
  vdt v_sum2a = vec_madd(v_xa, v_mhu_B2, v_mhu_B1);//BODYCOM0 
  vdt v_sum3a = vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);//BODYCOM0 
  vdt v_sum4a = vec_madd(v_xa, v_tauR_D2, v_tauR_D1);//BODYCOM0 
//BODYCOM0 
// 2  //BODYCOM0 
/*//BODYCOM0 
   sum1[ii] =  mhu_a[2]    + VM[ii]*sum1[ii];//BODYCOM0 
   sum1[ii] =  mhu_a[1]    + VM[ii]*sum1[ii];//BODYCOM0 
   sum2[ii] =  mhu_a[7]    + VM[ii]*sum2[ii];//BODYCOM0 
   sum2[ii] =  mhu_a[6]    + VM[ii]*sum2[ii];//BODYCOM0 
   sum3[ii] =  tauRdt_a[6] + VM[ii]*sum3[ii];//BODYCOM0 
   sum3[ii] =  tauRdt_a[5] + VM[ii]*sum3[ii];//BODYCOM0 
   sum4[ii] =  tauR_a[12]  + VM[ii]*sum4[ii]; //BODYCOM0 
   sum4[ii] =  tauR_a[11]  + VM[ii]*sum4[ii]; //BODYCOM0 
*///BODYCOM0 
  v_mhu_A1    = VEC_LDS(0, &mhu_a[1]);//BODYCOM0 
  v_mhu_A2    = VEC_LDS(0, &mhu_a[2]);//BODYCOM0 
  v_mhu_B1    = VEC_LDS(0, &mhu_a[6]);//BODYCOM0 
  v_mhu_B2    = VEC_LDS(0, &mhu_a[7]);//BODYCOM0 
  v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[5]);//BODYCOM0 
  v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[6]);//BODYCOM0 
  v_tauR_D1   = VEC_LDS(0, &tauR_a[11]);//BODYCOM0 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[12]);//BODYCOM0 
//BODYCOM0 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM0 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM0 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM0 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM0 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM0 
//BODYCOM0 
// 2  //BODYCOM0 
/*//BODYCOM0 
   sum1[ii] = mhu_a[0]    + VM[ii]*sum1[ii];//BODYCOM0 
   sum2[ii] = mhu_a[5]    + VM[ii]*sum2[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[4] + VM[ii]*sum3[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[3] + VM[ii]*sum3[ii];//BODYCOM0 
   sum4[ii] = tauR_a[10]  + VM[ii]*sum4[ii]; //BODYCOM0 
   sum4[ii] = tauR_a[9]   + VM[ii]*sum4[ii]; //BODYCOM0 
*///BODYCOM0 
//BODYCOM0 
   v_mhu_A1    = VEC_LDS(0, &mhu_a[0]);//BODYCOM0 
   v_mhu_B1    = VEC_LDS(0, &mhu_a[5]);//BODYCOM0 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[3]);//BODYCOM0 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[4]);//BODYCOM0 
   v_tauR_D1   = VEC_LDS(0, &tauR_a[9]);//BODYCOM0 
   v_tauR_D2   = VEC_LDS(0, &tauR_a[10]);//BODYCOM0 
//BODYCOM0 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM0 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM0 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM0 
//BODYCOM0 
// 3//BODYCOM0 
/*//BODYCOM0 
   sum3[ii] = tauRdt_a[2]  + VM[ii]*sum3[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[1]  + VM[ii]*sum3[ii];//BODYCOM0 
*///BODYCOM0 
//BODYCOM0 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[1]);//BODYCOM0 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[2]);//BODYCOM0 
//BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
 //BODYCOM0 
// 4//BODYCOM0 
/*//BODYCOM0 
   sum3[ii] = tauRdt_a[0]  + VM[ii]*sum3[ii];//BODYCOM0 
*///BODYCOM0 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[0]);//BODYCOM0 
 //BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
//BODYCOM0 
/*//BODYCOM0 
   double tauRdt= sum3/sum4; //BODYCOM0 
   double mhu= sum1/sum2;//BODYCOM0 
*///BODYCOM0 
   //BODYCOM0 
//BODYCOM0 
/*//BODYCOM0 
tauRdt[ii] = sum3[ii]/sum4[ii]; //BODYCOM0 
//BODYCOM0 
   mhu[ii]= sum1[ii]/sum2[ii];//BODYCOM0 
//BODYCOM0 
   g[ii] +=  mhu[ii]*tauRdt[ii] - g[ii]*tauRdt[ii];//BODYCOM0 
*///BODYCOM0 
// Still need to reuse the register I have above ... the code will look like gibberish because//BODYCOM0 
// the targets won't make sense.//BODYCOM0 
// A//BODYCOM0 
   v_xa      = vec_swdiv_nochk(v_sum3a, v_sum4a); //BODYCOM0 
   v_mhu_A2  = vec_swdiv_nochk(v_sum1a,v_sum2a); //BODYCOM0 
   v_mhu_B2  = vec_ld(0, &g[ii]);  //BODYCOM0 
   v_tauR_D2 = vec_sub(v_mhu_A2,   v_mhu_B2);  //BODYCOM0 
   v_mhu_B2  = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);//BODYCOM0 

   
   //vec_st(v_mhu_B2, 0, &g[ii]);
   switch (misalignment) {
   case 1:
     g[ii+0] = v_mhu_B2 get [1];
     g[ii+1] = v_mhu_B2 get [2];
     g[ii+2] = v_mhu_B2 get [3];
     break;

   case 2:
     g[ii+0] = v_mhu_B2 get [2];
     g[ii+1] = v_mhu_B2 get [3];
     break ;

   case 3:
     g[ii+0] = v_mhu_B2 get [3];
     break;

     default:
    //println("BAD MISALIGNMENT: %d\n", misalignment); fflush(stdout);
     assert(0);
   }

   ii = 4 - misalignment;
 }

 for (ii;ii<TTROUND(nCells,8);ii+=8)
 { 
   vdt v_xa   = vec_ld(0, &VM[ii]);  
   vdt v_xb   = vec_ld(0, &VM[ii+4]);  
 
// 1
 /*
   sum1[ii] =  mhu_a[4];
   sum1[ii] =  mhu_a[3]    + VM[ii]*sum1[ii];
   sum2[ii] =  mhu_a[9];
   sum2[ii] =  mhu_a[8]    + VM[ii]*sum2[ii];
   sum3[ii] =  tauRdt_a[8];
   sum3[ii] =  tauRdt_a[7] + VM[ii]*sum3[ii];
   sum4[ii] =  tauR_a[14]; 
   sum4[ii] =  tauR_a[13]  + VM[ii]*sum4[ii]; 
*/
  vdt v_mhu_A1    = VEC_LDS(0, &mhu_a[3]);
  vdt v_mhu_A2    = VEC_LDS(0, &mhu_a[4]);
  vdt v_mhu_B1    = VEC_LDS(0, &mhu_a[8]);
  vdt v_mhu_B2    = VEC_LDS(0, &mhu_a[9]);
  vdt v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[7]);
  vdt v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[8]);
  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[13]);
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[14]);

  vdt v_sum1a = vec_madd(v_xa, v_mhu_A2, v_mhu_A1);
  vdt v_sum1b = vec_madd(v_xb, v_mhu_A2, v_mhu_A1);
  vdt v_sum2a = vec_madd(v_xa, v_mhu_B2, v_mhu_B1);
  vdt v_sum2b = vec_madd(v_xb, v_mhu_B2, v_mhu_B1);
  vdt v_sum3a = vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);
  vdt v_sum3b = vec_madd(v_xb, v_tauRdt_C2, v_tauRdt_C1);
  vdt v_sum4a = vec_madd(v_xa, v_tauR_D2, v_tauR_D1);
  vdt v_sum4b = vec_madd(v_xb, v_tauR_D2, v_tauR_D1);

// 2  
/*
   sum1[ii] =  mhu_a[2]    + VM[ii]*sum1[ii];
   sum1[ii] =  mhu_a[1]    + VM[ii]*sum1[ii];
   sum2[ii] =  mhu_a[7]    + VM[ii]*sum2[ii];
   sum2[ii] =  mhu_a[6]    + VM[ii]*sum2[ii];
   sum3[ii] =  tauRdt_a[6] + VM[ii]*sum3[ii];
   sum3[ii] =  tauRdt_a[5] + VM[ii]*sum3[ii];
   sum4[ii] =  tauR_a[12]  + VM[ii]*sum4[ii]; 
   sum4[ii] =  tauR_a[11]  + VM[ii]*sum4[ii]; 
*/
  v_mhu_A1    = VEC_LDS(0, &mhu_a[1]);
  v_mhu_A2    = VEC_LDS(0, &mhu_a[2]);
  v_mhu_B1    = VEC_LDS(0, &mhu_a[6]);
  v_mhu_B2    = VEC_LDS(0, &mhu_a[7]);
  v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[5]);
  v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[6]);
  v_tauR_D1   = VEC_LDS(0, &tauR_a[11]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[12]);

   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);

   v_sum1b = vec_madd(v_xb, v_sum1b, v_mhu_A2);
   v_sum1b = vec_madd(v_xb, v_sum1b, v_mhu_A1);
   v_sum2b = vec_madd(v_xb, v_sum2b, v_mhu_B2);
   v_sum2b = vec_madd(v_xb, v_sum2b, v_mhu_B1); 
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   v_sum4b = vec_madd(v_xb, v_sum4b, v_tauR_D2);
   v_sum4b = vec_madd(v_xb, v_sum4b, v_tauR_D1);


// 2  
/*
   sum1[ii] = mhu_a[0]    + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[5]    + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[4] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[3] + VM[ii]*sum3[ii];
   sum4[ii] = tauR_a[10]  + VM[ii]*sum4[ii]; 
   sum4[ii] = tauR_a[9]   + VM[ii]*sum4[ii]; 
*/

   v_mhu_A1    = VEC_LDS(0, &mhu_a[0]);
   v_mhu_B1    = VEC_LDS(0, &mhu_a[5]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[3]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[4]);
   v_tauR_D1   = VEC_LDS(0, &tauR_a[9]);
   v_tauR_D2   = VEC_LDS(0, &tauR_a[10]);

   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);

   v_sum1b = vec_madd(v_xb, v_sum1b, v_mhu_A1);
   v_sum2b = vec_madd(v_xb, v_sum2b, v_mhu_B1); 
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   v_sum4b = vec_madd(v_xb, v_sum4b, v_tauR_D2);
   v_sum4b = vec_madd(v_xb, v_sum4b, v_tauR_D1);


// 3
/*
   sum3[ii] = tauRdt_a[2]  + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[1]  + VM[ii]*sum3[ii];
*/

   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[1]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[2]);

   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
 
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
 

// 4
/*
   sum3[ii] = tauRdt_a[0]  + VM[ii]*sum3[ii];
*/
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[0]);
 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C1);

/*
   double tauRdt= sum3/sum4; 
   double mhu= sum1/sum2;
*/
   

/*
tauRdt[ii] = sum3[ii]/sum4[ii]; 

   mhu[ii]= sum1[ii]/sum2[ii];

   g[ii] +=  mhu[ii]*tauRdt[ii] - g[ii]*tauRdt[ii];
*/
// Still need to reuse the register I have above ... the code will look like gibberish because
// the targets won't make sense.
// A
   v_xa      = vec_swdiv_nochk(v_sum3a, v_sum4a); 
   v_mhu_A2  = vec_swdiv_nochk(v_sum1a,v_sum2a); 
   v_mhu_B2  = vec_ld(0, &g[ii]);  
   v_tauR_D2 = vec_sub(v_mhu_A2,   v_mhu_B2);  
   v_mhu_B2  = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);
   vec_st(v_mhu_B2, 0, &g[ii]);

// B
   v_xb      = vec_swdiv_nochk(v_sum3b, v_sum4b); 
   v_mhu_A1  = vec_swdiv_nochk(v_sum1b, v_sum2b); 
   v_mhu_B1  = vec_ld(0, &g[ii+4]);  
   v_tauR_D1 = vec_sub(v_mhu_A1, v_mhu_B1);  
   v_mhu_B1  = vec_madd(v_tauR_D1, v_xb, v_mhu_B1);
   vec_st(v_mhu_B1, 0, &g[ii+4]);

  }
 for (ii;ii<TTROUND(nCells,4);ii+=4)
 { 
 //BODYSTART
   vdt v_xa   = vec_ld(0, &VM[ii]);  
 
// 1
 /*
   sum1[ii] =  mhu_a[4];
   sum1[ii] =  mhu_a[3]    + VM[ii]*sum1[ii];
   sum2[ii] =  mhu_a[9];
   sum2[ii] =  mhu_a[8]    + VM[ii]*sum2[ii];
   sum3[ii] =  tauRdt_a[8];
   sum3[ii] =  tauRdt_a[7] + VM[ii]*sum3[ii];
   sum4[ii] =  tauR_a[14]; 
   sum4[ii] =  tauR_a[13]  + VM[ii]*sum4[ii]; 
*/
  vdt v_mhu_A1    = VEC_LDS(0, &mhu_a[3]);
  vdt v_mhu_A2    = VEC_LDS(0, &mhu_a[4]);
  vdt v_mhu_B1    = VEC_LDS(0, &mhu_a[8]);
  vdt v_mhu_B2    = VEC_LDS(0, &mhu_a[9]);
  vdt v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[7]);
  vdt v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[8]);
  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[13]);
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[14]);

  vdt v_sum1a = vec_madd(v_xa, v_mhu_A2, v_mhu_A1);
  vdt v_sum2a = vec_madd(v_xa, v_mhu_B2, v_mhu_B1);
  vdt v_sum3a = vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);
  vdt v_sum4a = vec_madd(v_xa, v_tauR_D2, v_tauR_D1);

// 2  
/*
   sum1[ii] =  mhu_a[2]    + VM[ii]*sum1[ii];
   sum1[ii] =  mhu_a[1]    + VM[ii]*sum1[ii];
   sum2[ii] =  mhu_a[7]    + VM[ii]*sum2[ii];
   sum2[ii] =  mhu_a[6]    + VM[ii]*sum2[ii];
   sum3[ii] =  tauRdt_a[6] + VM[ii]*sum3[ii];
   sum3[ii] =  tauRdt_a[5] + VM[ii]*sum3[ii];
   sum4[ii] =  tauR_a[12]  + VM[ii]*sum4[ii]; 
   sum4[ii] =  tauR_a[11]  + VM[ii]*sum4[ii]; 
*/
  v_mhu_A1    = VEC_LDS(0, &mhu_a[1]);
  v_mhu_A2    = VEC_LDS(0, &mhu_a[2]);
  v_mhu_B1    = VEC_LDS(0, &mhu_a[6]);
  v_mhu_B2    = VEC_LDS(0, &mhu_a[7]);
  v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[5]);
  v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[6]);
  v_tauR_D1   = VEC_LDS(0, &tauR_a[11]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[12]);

   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);

// 2  
/*
   sum1[ii] = mhu_a[0]    + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[5]    + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[4] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[3] + VM[ii]*sum3[ii];
   sum4[ii] = tauR_a[10]  + VM[ii]*sum4[ii]; 
   sum4[ii] = tauR_a[9]   + VM[ii]*sum4[ii]; 
*/

   v_mhu_A1    = VEC_LDS(0, &mhu_a[0]);
   v_mhu_B1    = VEC_LDS(0, &mhu_a[5]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[3]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[4]);
   v_tauR_D1   = VEC_LDS(0, &tauR_a[9]);
   v_tauR_D2   = VEC_LDS(0, &tauR_a[10]);

   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);

// 3
/*
   sum3[ii] = tauRdt_a[2]  + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[1]  + VM[ii]*sum3[ii];
*/

   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[1]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[2]);

   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
 
// 4
/*
   sum3[ii] = tauRdt_a[0]  + VM[ii]*sum3[ii];
*/
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[0]);
 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);

/*
   double tauRdt= sum3/sum4; 
   double mhu= sum1/sum2;
*/
   

/*
tauRdt[ii] = sum3[ii]/sum4[ii]; 

   mhu[ii]= sum1[ii]/sum2[ii];

   g[ii] +=  mhu[ii]*tauRdt[ii] - g[ii]*tauRdt[ii];
*/
// Still need to reuse the register I have above ... the code will look like gibberish because
// the targets won't make sense.
// A
   v_xa      = vec_swdiv_nochk(v_sum3a, v_sum4a); 
   v_mhu_A2  = vec_swdiv_nochk(v_sum1a,v_sum2a); 
   v_mhu_B2  = vec_ld(0, &g[ii]);  
   v_tauR_D2 = vec_sub(v_mhu_A2,   v_mhu_B2);  
   v_mhu_B2  = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);
   vec_st(v_mhu_B2, 0, &g[ii]); //BODYEND
  }

 // the elements in the end of the vector should be ignored
 int leftover = nCells - ii;
 if (leftover) { 
      vdt v_xa   = vec_ld(0, &VM[ii]);  //BODYCOM1 
 //BODYCOM1 
// 1//BODYCOM1 
 /*//BODYCOM1 
   sum1[ii] =  mhu_a[4];//BODYCOM1 
   sum1[ii] =  mhu_a[3]    + VM[ii]*sum1[ii];//BODYCOM1 
   sum2[ii] =  mhu_a[9];//BODYCOM1 
   sum2[ii] =  mhu_a[8]    + VM[ii]*sum2[ii];//BODYCOM1 
   sum3[ii] =  tauRdt_a[8];//BODYCOM1 
   sum3[ii] =  tauRdt_a[7] + VM[ii]*sum3[ii];//BODYCOM1 
   sum4[ii] =  tauR_a[14]; //BODYCOM1 
   sum4[ii] =  tauR_a[13]  + VM[ii]*sum4[ii]; //BODYCOM1 
*///BODYCOM1 
  vdt v_mhu_A1    = VEC_LDS(0, &mhu_a[3]);//BODYCOM1 
  vdt v_mhu_A2    = VEC_LDS(0, &mhu_a[4]);//BODYCOM1 
  vdt v_mhu_B1    = VEC_LDS(0, &mhu_a[8]);//BODYCOM1 
  vdt v_mhu_B2    = VEC_LDS(0, &mhu_a[9]);//BODYCOM1 
  vdt v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[7]);//BODYCOM1 
  vdt v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[8]);//BODYCOM1 
  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[13]);//BODYCOM1 
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[14]);//BODYCOM1 
//BODYCOM1 
  vdt v_sum1a = vec_madd(v_xa, v_mhu_A2, v_mhu_A1);//BODYCOM1 
  vdt v_sum2a = vec_madd(v_xa, v_mhu_B2, v_mhu_B1);//BODYCOM1 
  vdt v_sum3a = vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);//BODYCOM1 
  vdt v_sum4a = vec_madd(v_xa, v_tauR_D2, v_tauR_D1);//BODYCOM1 
//BODYCOM1 
// 2  //BODYCOM1 
/*//BODYCOM1 
   sum1[ii] =  mhu_a[2]    + VM[ii]*sum1[ii];//BODYCOM1 
   sum1[ii] =  mhu_a[1]    + VM[ii]*sum1[ii];//BODYCOM1 
   sum2[ii] =  mhu_a[7]    + VM[ii]*sum2[ii];//BODYCOM1 
   sum2[ii] =  mhu_a[6]    + VM[ii]*sum2[ii];//BODYCOM1 
   sum3[ii] =  tauRdt_a[6] + VM[ii]*sum3[ii];//BODYCOM1 
   sum3[ii] =  tauRdt_a[5] + VM[ii]*sum3[ii];//BODYCOM1 
   sum4[ii] =  tauR_a[12]  + VM[ii]*sum4[ii]; //BODYCOM1 
   sum4[ii] =  tauR_a[11]  + VM[ii]*sum4[ii]; //BODYCOM1 
*///BODYCOM1 
  v_mhu_A1    = VEC_LDS(0, &mhu_a[1]);//BODYCOM1 
  v_mhu_A2    = VEC_LDS(0, &mhu_a[2]);//BODYCOM1 
  v_mhu_B1    = VEC_LDS(0, &mhu_a[6]);//BODYCOM1 
  v_mhu_B2    = VEC_LDS(0, &mhu_a[7]);//BODYCOM1 
  v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[5]);//BODYCOM1 
  v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[6]);//BODYCOM1 
  v_tauR_D1   = VEC_LDS(0, &tauR_a[11]);//BODYCOM1 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[12]);//BODYCOM1 
//BODYCOM1 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM1 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM1 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM1 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM1 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM1 
//BODYCOM1 
// 2  //BODYCOM1 
/*//BODYCOM1 
   sum1[ii] = mhu_a[0]    + VM[ii]*sum1[ii];//BODYCOM1 
   sum2[ii] = mhu_a[5]    + VM[ii]*sum2[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[4] + VM[ii]*sum3[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[3] + VM[ii]*sum3[ii];//BODYCOM1 
   sum4[ii] = tauR_a[10]  + VM[ii]*sum4[ii]; //BODYCOM1 
   sum4[ii] = tauR_a[9]   + VM[ii]*sum4[ii]; //BODYCOM1 
*///BODYCOM1 
//BODYCOM1 
   v_mhu_A1    = VEC_LDS(0, &mhu_a[0]);//BODYCOM1 
   v_mhu_B1    = VEC_LDS(0, &mhu_a[5]);//BODYCOM1 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[3]);//BODYCOM1 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[4]);//BODYCOM1 
   v_tauR_D1   = VEC_LDS(0, &tauR_a[9]);//BODYCOM1 
   v_tauR_D2   = VEC_LDS(0, &tauR_a[10]);//BODYCOM1 
//BODYCOM1 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM1 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM1 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM1 
//BODYCOM1 
// 3//BODYCOM1 
/*//BODYCOM1 
   sum3[ii] = tauRdt_a[2]  + VM[ii]*sum3[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[1]  + VM[ii]*sum3[ii];//BODYCOM1 
*///BODYCOM1 
//BODYCOM1 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[1]);//BODYCOM1 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[2]);//BODYCOM1 
//BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
 //BODYCOM1 
// 4//BODYCOM1 
/*//BODYCOM1 
   sum3[ii] = tauRdt_a[0]  + VM[ii]*sum3[ii];//BODYCOM1 
*///BODYCOM1 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[0]);//BODYCOM1 
 //BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
//BODYCOM1 
/*//BODYCOM1 
   double tauRdt= sum3/sum4; //BODYCOM1 
   double mhu= sum1/sum2;//BODYCOM1 
*///BODYCOM1 
   //BODYCOM1 
//BODYCOM1 
/*//BODYCOM1 
tauRdt[ii] = sum3[ii]/sum4[ii]; //BODYCOM1 
//BODYCOM1 
   mhu[ii]= sum1[ii]/sum2[ii];//BODYCOM1 
//BODYCOM1 
   g[ii] +=  mhu[ii]*tauRdt[ii] - g[ii]*tauRdt[ii];//BODYCOM1 
*///BODYCOM1 
// Still need to reuse the register I have above ... the code will look like gibberish because//BODYCOM1 
// the targets won't make sense.//BODYCOM1 
// A//BODYCOM1 
   v_xa      = vec_swdiv_nochk(v_sum3a, v_sum4a); //BODYCOM1 
   v_mhu_A2  = vec_swdiv_nochk(v_sum1a,v_sum2a); //BODYCOM1 
   v_mhu_B2  = vec_ld(0, &g[ii]);  //BODYCOM1 
   v_tauR_D2 = vec_sub(v_mhu_A2,   v_mhu_B2);  //BODYCOM1 
   v_mhu_B2  = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);//BODYCOM1 

   
   //vec_st(v_mhu_B2, 0, &g[ii]);
   switch (leftover) {
   case 1:
     g[ii+0] = v_mhu_B2 get [0];
     break;

   case 2:
     g[ii+0] = v_mhu_B2 get [0];
     g[ii+1] = v_mhu_B2 get [1];
     break;

   case 3:
     g[ii+0] = v_mhu_B2 get [0];
     g[ii+1] = v_mhu_B2 get [1];
     g[ii+2] = v_mhu_B2 get [2];
     break;

     default:
    //println("BAD LEFTOVER: %d\n", leftover); fflush(stdout);
    assert(0);
   }
 }

}




// 4134
void update_rGate_v1(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a)
{
  //update_rGate(dt, nCells, VM, g, mhu_a, tauR_a); return;
  typedef vector4double vdt;

int gateIndex=6;
int mhu_l=5 ;
int mhu_m= 8; // Was 5
int tauR_l= 5;
int tauR_m= 7;
   
 double  tauRdt_a[7];


         tauRdt_a[0]  = tauR_a[0]*dt;
         tauRdt_a[1]  = tauR_a[1]*dt;
         tauRdt_a[2]  = tauR_a[2]*dt;
         tauRdt_a[3]  = tauR_a[3]*dt;
         tauRdt_a[4]  = tauR_a[4]*dt;
         tauRdt_a[5]  = tauR_a[5]*dt;
         tauRdt_a[6]  = tauR_a[6]*dt; 
 int ii = 0;

 // 0 aligned
 // 1: load produces (JUNK, &VM[0], &VM[1], &VM[2])
 // 2: load produces (JUNK, JUNK, &VM[0], &VM[1])
 // 3: load produces (JUNK, JUNK, JUNK, &VM[0])
 int misalignment = ((long int)&VM[0] % 32) / 8;
 assert(misalignment == 0);

 if (misalignment) { 
      vdt v_xa   = vec_ld(0, &VM[ii]);  //BODYCOM0 
//BODYCOM0 
// 1//BODYCOM0 
  vdt v_mhu_A1  = VEC_LDS(0, &mhu_a[6]);//BODYCOM0 
  vdt v_mhu_A2  = VEC_LDS(0, &mhu_a[7]);//BODYCOM0 
  vdt v_mhu_B1 = VEC_LDS(0, &mhu_a[11]);//BODYCOM0 
  vdt v_mhu_B2 = VEC_LDS(0, &mhu_a[12]);//BODYCOM0 
  vdt v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[5]);//BODYCOM0 
  vdt v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[6]);//BODYCOM0 
  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[10]);//BODYCOM0 
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[11]);//BODYCOM0 
//BODYCOM0 
   vdt v_sum1a =  vec_madd(v_xa, v_mhu_A2, v_mhu_A1);//BODYCOM0 
   vdt v_sum2a =  vec_madd(v_xa, v_mhu_B2, v_mhu_B1);//BODYCOM0 
   vdt v_sum3a =  vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);//BODYCOM0 
   vdt v_sum4a =  vec_madd(v_xa, v_tauR_D2, v_tauR_D1);//BODYCOM0 
// 2//BODYCOM0 
   v_mhu_A1  = VEC_LDS(0, &mhu_a[4]);//BODYCOM0 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[5]);//BODYCOM0 
   v_mhu_B1 = VEC_LDS(0, &mhu_a[9]);//BODYCOM0 
   v_mhu_B2 = VEC_LDS(0, &mhu_a[10]);//BODYCOM0 
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[3]);//BODYCOM0 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[4]);//BODYCOM0 
  v_tauR_D1   = VEC_LDS(0, &tauR_a[8]);//BODYCOM0 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[9]);//BODYCOM0 
//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM0 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM0 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM0 
// 3//BODYCOM0 
//BODYCOM0 
   v_mhu_A1  = VEC_LDS(0, &mhu_a[2]);//BODYCOM0 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[3]);//BODYCOM0 
   v_mhu_B2 = VEC_LDS(0, &mhu_a[8]);//BODYCOM0 
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[1]);//BODYCOM0 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[2]);//BODYCOM0 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[7]);//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM0 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM0 
// 4//BODYCOM0 
//BODYCOM0 
   v_mhu_A1  = VEC_LDS(0, &mhu_a[0]);//BODYCOM0 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[1]);//BODYCOM0 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[0]);//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
//BODYCOM0 
// Still need to reuse the register I have above ... the code will look like gibberish because//BODYCOM0 
// the targets won't make sense.//BODYCOM0 
// A//BODYCOM0 
   vdt temp1a   = vec_swdiv_nochk(v_sum3a,v_sum4a); //BODYCOM0 
   vdt temp2a   = vec_swdiv_nochk(v_sum1a,v_sum2a); //BODYCOM0 
   vdt temp3a   = vec_ld(0, &g[ii]);  //BODYCOM0 
   temp2a   = vec_sub(temp2a, temp3a);  //BODYCOM0 
   temp3a   = vec_madd(temp2a, temp1a, temp3a);//BODYCOM0 

   
   //vec_st(temp3a, 0, &g[ii]);
   switch (misalignment) {
   case 1:
     g[ii+0] = temp3a get [1];
     g[ii+1] = temp3a get [2];
     g[ii+2] = temp3a get [3];
     break;

   case 2:
     g[ii+0] = temp3a get [2];
     g[ii+1] = temp3a get [3];
     break;

   case 3:
     g[ii+0] = temp3a get [3];
     break;

     default:
    //println("BAD MISALIGNMENT: %d\n", misalignment); fflush(stdout);
     assert(0);
   }

   ii = 4 - misalignment;
 }

 for (ii;ii<TTROUND(nCells,16);ii+=16)
 { 
   vdt v_xa   = vec_ld(0, &VM[ii]);  
   vdt v_xb   = vec_ld(0, &VM[ii+4]);  
   vdt v_xc   = vec_ld(0, &VM[ii+8]);  
   vdt v_xd   = vec_ld(0, &VM[ii+12]);  

// 1
  vdt v_mhu_A1  = VEC_LDS(0, &mhu_a[6]);
  vdt v_mhu_A2  = VEC_LDS(0, &mhu_a[7]);
  vdt v_mhu_B1 = VEC_LDS(0, &mhu_a[11]);
  vdt v_mhu_B2 = VEC_LDS(0, &mhu_a[12]);
  vdt v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[5]);
  vdt v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[6]);
  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[10]);
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[11]);

   vdt v_sum1a =  vec_madd(v_xa, v_mhu_A2, v_mhu_A1);
   vdt v_sum1b =  vec_madd(v_xb, v_mhu_A2, v_mhu_A1);
   vdt v_sum1c =  vec_madd(v_xc, v_mhu_A2, v_mhu_A1);
   vdt v_sum1d =  vec_madd(v_xd, v_mhu_A2, v_mhu_A1);

   vdt v_sum2a =  vec_madd(v_xa, v_mhu_B2, v_mhu_B1);
   vdt v_sum2b =  vec_madd(v_xb, v_mhu_B2, v_mhu_B1);
   vdt v_sum2c =  vec_madd(v_xc, v_mhu_B2, v_mhu_B1);
   vdt v_sum2d =  vec_madd(v_xd, v_mhu_B2, v_mhu_B1);

   vdt v_sum3a =  vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);
   vdt v_sum3b =  vec_madd(v_xb, v_tauRdt_C2, v_tauRdt_C1);
   vdt v_sum3c =  vec_madd(v_xc, v_tauRdt_C2, v_tauRdt_C1);
   vdt v_sum3d =  vec_madd(v_xd, v_tauRdt_C2, v_tauRdt_C1);

   vdt v_sum4a =  vec_madd(v_xa, v_tauR_D2, v_tauR_D1);
   vdt v_sum4b =  vec_madd(v_xb, v_tauR_D2, v_tauR_D1);
   vdt v_sum4c =  vec_madd(v_xc, v_tauR_D2, v_tauR_D1);
   vdt v_sum4d =  vec_madd(v_xd, v_tauR_D2, v_tauR_D1);
// 2
   v_mhu_A1  = VEC_LDS(0, &mhu_a[4]);
   v_mhu_A2  = VEC_LDS(0, &mhu_a[5]);
   v_mhu_B1 = VEC_LDS(0, &mhu_a[9]);
   v_mhu_B2 = VEC_LDS(0, &mhu_a[10]);
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[3]);
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[4]);
  v_tauR_D1   = VEC_LDS(0, &tauR_a[8]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[9]);

   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A2);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A1);
   v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B2);
   v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B1); 
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D2);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D1);
   v_sum1c =  vec_madd(v_xc, v_sum1c, v_mhu_A2);
   v_sum1c =  vec_madd(v_xc, v_sum1c, v_mhu_A1);
   v_sum2c =  vec_madd(v_xc, v_sum2c, v_mhu_B2);
   v_sum2c =  vec_madd(v_xc, v_sum2c, v_mhu_B1); 
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C2);
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C1);
   v_sum4c =  vec_madd(v_xc, v_sum4c, v_tauR_D2);
   v_sum4c =  vec_madd(v_xc, v_sum4c, v_tauR_D1);
   v_sum1d =  vec_madd(v_xd, v_sum1d, v_mhu_A2);
   v_sum1d =  vec_madd(v_xd, v_sum1d, v_mhu_A1);
   v_sum2d =  vec_madd(v_xd, v_sum2d, v_mhu_B2);
   v_sum2d =  vec_madd(v_xd, v_sum2d, v_mhu_B1); 
   v_sum3d =  vec_madd(v_xd, v_sum3d, v_tauRdt_C2);
   v_sum3d =  vec_madd(v_xd, v_sum3d, v_tauRdt_C1);
   v_sum4d =  vec_madd(v_xd, v_sum4d, v_tauR_D2);
   v_sum4d =  vec_madd(v_xd, v_sum4d, v_tauR_D1);
// 3

   v_mhu_A1  = VEC_LDS(0, &mhu_a[2]);
   v_mhu_A2  = VEC_LDS(0, &mhu_a[3]);
   v_mhu_B2 = VEC_LDS(0, &mhu_a[8]);
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[1]);
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[2]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[7]);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A2);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A1);
   v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B2);
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D2);
   v_sum1c =  vec_madd(v_xc, v_sum1c, v_mhu_A2);
   v_sum1c =  vec_madd(v_xc, v_sum1c, v_mhu_A1);
   v_sum2c =  vec_madd(v_xc, v_sum2c, v_mhu_B2);
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C2);
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C1);
   v_sum4c =  vec_madd(v_xc, v_sum4c, v_tauR_D2);
   v_sum1d =  vec_madd(v_xd, v_sum1d, v_mhu_A2);
   v_sum1d =  vec_madd(v_xd, v_sum1d, v_mhu_A1);
   v_sum2d =  vec_madd(v_xd, v_sum2d, v_mhu_B2);
   v_sum3d =  vec_madd(v_xd, v_sum3d, v_tauRdt_C2);
   v_sum3d =  vec_madd(v_xd, v_sum3d, v_tauRdt_C1);
   v_sum4d =  vec_madd(v_xd, v_sum4d, v_tauR_D2);
// 4

   v_mhu_A1  = VEC_LDS(0, &mhu_a[0]);
   v_mhu_A2  = VEC_LDS(0, &mhu_a[1]);
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[0]);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A2);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A1);
   v_sum1c =  vec_madd(v_xc, v_sum1c, v_mhu_A2);
   v_sum1c =  vec_madd(v_xc, v_sum1c, v_mhu_A1);
   v_sum1d =  vec_madd(v_xd, v_sum1d, v_mhu_A2);
   v_sum1d =  vec_madd(v_xd, v_sum1d, v_mhu_A1);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C2);
   v_sum3d =  vec_madd(v_xd, v_sum3d, v_tauRdt_C2);

// Still need to reuse the register I have above ... the code will look like gibberish because
// the targets won't make sense.
// A
   vdt temp1a   = vec_swdiv_nochk(v_sum3a,v_sum4a); 
   vdt temp2a   = vec_swdiv_nochk(v_sum1a,v_sum2a); 
   vdt temp3a   = vec_ld(0, &g[ii]);  
   temp2a   = vec_sub(temp2a, temp3a);  
   temp3a   = vec_madd(temp2a, temp1a, temp3a);
   vec_st(temp3a, 0, &g[ii]);
   /* This can also be implemented as the following, which has a higher flop rate and is more straightforward
   temp1a  = vec_swdiv_nochk(v_sum3a,v_sum4a); 
   temp2a  = vec_swdiv_nochk(v_sum1a,v_sum2a); 
   temp3a  = vec_ld(0, &g[ii]);
   temp2a   = vec_madd(temp1a, temp2a, temp3a);
   temp3a  = vec_nmadd(temp1a, temp3a, temp2a );
   vec_st(temp3a, 0, &g[ii]); */
// B
   vdt temp1b   = vec_swdiv_nochk(v_sum3b,v_sum4b); 
   vdt temp2b   = vec_swdiv_nochk(v_sum1b,v_sum2b); 
   vdt temp3b   = vec_ld(0, &g[ii+4]);  
   temp2b   = vec_sub(temp2b, temp3b);  
   temp3b   = vec_madd(temp2b, temp1b, temp3b);
   vec_st(temp3b, 0, &g[ii+4]);
// C
   vdt temp1c   = vec_swdiv_nochk(v_sum3c,v_sum4c); 
   vdt temp2c   = vec_swdiv_nochk(v_sum1c,v_sum2c); 
   vdt temp3c   = vec_ld(0, &g[ii+8]);  
   temp2c   = vec_sub(temp2c, temp3c);  
   temp3c   = vec_madd(temp2c, temp1c, temp3c);
   vec_st(temp3c, 0, &g[ii+8]);
// D
   vdt temp1d   = vec_swdiv_nochk(v_sum3d,v_sum4d); 
   vdt temp2d   = vec_swdiv_nochk(v_sum1d,v_sum2d); 
   vdt temp3d   = vec_ld(0, &g[ii+12]);  
   temp2d   = vec_sub(temp2d, temp3d);  
   temp3d   = vec_madd(temp2d, temp1d, temp3d);
   vec_st(temp3d, 0, &g[ii+12]);
 }
 for (ii;ii<TTROUND(nCells,4);ii+=4)
 { 
 //BODYSTART
   vdt v_xa   = vec_ld(0, &VM[ii]);  

// 1
  vdt v_mhu_A1  = VEC_LDS(0, &mhu_a[6]);
  vdt v_mhu_A2  = VEC_LDS(0, &mhu_a[7]);
  vdt v_mhu_B1 = VEC_LDS(0, &mhu_a[11]);
  vdt v_mhu_B2 = VEC_LDS(0, &mhu_a[12]);
  vdt v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[5]);
  vdt v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[6]);
  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[10]);
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[11]);

   vdt v_sum1a =  vec_madd(v_xa, v_mhu_A2, v_mhu_A1);
   vdt v_sum2a =  vec_madd(v_xa, v_mhu_B2, v_mhu_B1);
   vdt v_sum3a =  vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);
   vdt v_sum4a =  vec_madd(v_xa, v_tauR_D2, v_tauR_D1);
// 2
   v_mhu_A1  = VEC_LDS(0, &mhu_a[4]);
   v_mhu_A2  = VEC_LDS(0, &mhu_a[5]);
   v_mhu_B1 = VEC_LDS(0, &mhu_a[9]);
   v_mhu_B2 = VEC_LDS(0, &mhu_a[10]);
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[3]);
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[4]);
  v_tauR_D1   = VEC_LDS(0, &tauR_a[8]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[9]);

   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);
// 3

   v_mhu_A1  = VEC_LDS(0, &mhu_a[2]);
   v_mhu_A2  = VEC_LDS(0, &mhu_a[3]);
   v_mhu_B2 = VEC_LDS(0, &mhu_a[8]);
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[1]);
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[2]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[7]);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   //v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
// 4

   v_mhu_A1  = VEC_LDS(0, &mhu_a[0]);
   v_mhu_A2  = VEC_LDS(0, &mhu_a[1]);
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[0]);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);

// Still need to reuse the register I have above ... the code will look like gibberish because
// the targets won't make sense.
// A
   vdt temp1a   = vec_swdiv_nochk(v_sum3a,v_sum4a); 
   vdt temp2a   = vec_swdiv_nochk(v_sum1a,v_sum2a); 
   vdt temp3a   = vec_ld(0, &g[ii]);  
   temp2a   = vec_sub(temp2a, temp3a);  
   temp3a   = vec_madd(temp2a, temp1a, temp3a);
   vec_st(temp3a, 0, &g[ii]); //BODYEND
   /* This can also be implemented as the following, which has a higher flop rate and is more straightforward
   temp1a  = vec_swdiv_nochk(v_sum3a,v_sum4a); 
   temp2a  = vec_swdiv_nochk(v_sum1a,v_sum2a); 
   temp3a  = vec_ld(0, &g[ii]);
   temp2a   = vec_madd(temp1a, temp2a, temp3a);
   temp3a  = vec_nmadd(temp1a, temp3a, temp2a );
   vec_st(temp3a, 0, &g[ii]); */
 }

 // the elements in the end of the vector should be ignored
 int leftover = nCells - ii;
 if (leftover) { 
      vdt v_xa   = vec_ld(0, &VM[ii]);  //BODYCOM1 
//BODYCOM1 
// 1//BODYCOM1 
  vdt v_mhu_A1  = VEC_LDS(0, &mhu_a[6]);//BODYCOM1 
  vdt v_mhu_A2  = VEC_LDS(0, &mhu_a[7]);//BODYCOM1 
  vdt v_mhu_B1 = VEC_LDS(0, &mhu_a[11]);//BODYCOM1 
  vdt v_mhu_B2 = VEC_LDS(0, &mhu_a[12]);//BODYCOM1 
  vdt v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[5]);//BODYCOM1 
  vdt v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[6]);//BODYCOM1 
  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[10]);//BODYCOM1 
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[11]);//BODYCOM1 
//BODYCOM1 
   vdt v_sum1a =  vec_madd(v_xa, v_mhu_A2, v_mhu_A1);//BODYCOM1 
   vdt v_sum2a =  vec_madd(v_xa, v_mhu_B2, v_mhu_B1);//BODYCOM1 
   vdt v_sum3a =  vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);//BODYCOM1 
   vdt v_sum4a =  vec_madd(v_xa, v_tauR_D2, v_tauR_D1);//BODYCOM1 
// 2//BODYCOM1 
   v_mhu_A1  = VEC_LDS(0, &mhu_a[4]);//BODYCOM1 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[5]);//BODYCOM1 
   v_mhu_B1 = VEC_LDS(0, &mhu_a[9]);//BODYCOM1 
   v_mhu_B2 = VEC_LDS(0, &mhu_a[10]);//BODYCOM1 
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[3]);//BODYCOM1 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[4]);//BODYCOM1 
  v_tauR_D1   = VEC_LDS(0, &tauR_a[8]);//BODYCOM1 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[9]);//BODYCOM1 
//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM1 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM1 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM1 
// 3//BODYCOM1 
//BODYCOM1 
   v_mhu_A1  = VEC_LDS(0, &mhu_a[2]);//BODYCOM1 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[3]);//BODYCOM1 
   v_mhu_B2 = VEC_LDS(0, &mhu_a[8]);//BODYCOM1 
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[1]);//BODYCOM1 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[2]);//BODYCOM1 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[7]);//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM1 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM1 
// 4//BODYCOM1 
//BODYCOM1 
   v_mhu_A1  = VEC_LDS(0, &mhu_a[0]);//BODYCOM1 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[1]);//BODYCOM1 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[0]);//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
//BODYCOM1 
// Still need to reuse the register I have above ... the code will look like gibberish because//BODYCOM1 
// the targets won't make sense.//BODYCOM1 
// A//BODYCOM1 
   vdt temp1a   = vec_swdiv_nochk(v_sum3a,v_sum4a); //BODYCOM1 
   vdt temp2a   = vec_swdiv_nochk(v_sum1a,v_sum2a); //BODYCOM1 
   vdt temp3a   = vec_ld(0, &g[ii]);  //BODYCOM1 
   temp2a   = vec_sub(temp2a, temp3a);  //BODYCOM1 
   temp3a   = vec_madd(temp2a, temp1a, temp3a);//BODYCOM1 

   
   //vec_st(temp3a, 0, &g[ii]);
   switch (leftover) {
   case 1:
     g[ii+0] = temp3a get [0];
     break;

   case 2:
     g[ii+0] = temp3a get [0];
     g[ii+1] = temp3a get [1];
     break;

   case 3:
     g[ii+0] = temp3a get [0];
     g[ii+1] = temp3a get [1];
     g[ii+2] = temp3a get [2];
     break;

     default:
    //println("BAD LEFTOVER: %d\n", leftover); fflush(stdout);
    assert(0);
   }
 }

}



void update_dGate_v1(double dt, int nCells,  double *VM, double *g, double *mhu_a, double *tauR_a)
{
  //update_dGate(dt, nCells,  VM, g, mhu_a, tauR_a); return;
  typedef vector4double vdt;

int gateIndex= 7;
int mhu_l= 5 ;
int mhu_m= 7;
int tauR_l=7 ;
int tauR_m=10;

/*
 int mhu_k  = 11;  
 int tauR_k = 16;  
*/
   
 double  tauRdt_a[10];

//  for (int j=9;j>=0;j--)    tauRdt_a[j] = tauR_a[j]*dt;
         tauRdt_a[0]  = tauR_a[0]*dt;
         tauRdt_a[1]  = tauR_a[1]*dt;
         tauRdt_a[2]  = tauR_a[2]*dt;
         tauRdt_a[3]  = tauR_a[3]*dt;
         tauRdt_a[4]  = tauR_a[4]*dt;
         tauRdt_a[5]  = tauR_a[5]*dt;
         tauRdt_a[6]  = tauR_a[6]*dt; 
         tauRdt_a[7]  = tauR_a[7]*dt;
         tauRdt_a[8]  = tauR_a[8]*dt;
         tauRdt_a[9]  = tauR_a[9]*dt;

 vdt v_mhu_a5  = VEC_LDS(0, &mhu_a[5]);
 vdt v_mhu_a6  = VEC_LDS(0, &mhu_a[6]); 
 vdt v_mhu_a10 = VEC_LDS(0, &mhu_a[10]);
 vdt v_mhu_a11 = VEC_LDS(0, &mhu_a[11]);   

 int ii=0;

 // 0 aligned
 // 1: load produces (JUNK, &VM[0], &VM[1], &VM[2])
 // 2: load produces (JUNK, JUNK, &VM[0], &VM[1])
 // 3: load produces (JUNK, JUNK, JUNK, &VM[0])
 int misalignment = ((long int)&VM[0] % 32) / 8;
 assert(misalignment == 0);

 if (misalignment) { 
      vdt v_xa   = vec_ld(0, &VM[ii]);  //BODYCOM0 
//BODYCOM0 
// 1//BODYCOM0 
/*//BODYCOM0 
   sum1[ii] = mhu_a[6];//BODYCOM0 
   sum1[ii] = mhu_a[5]    + VM[ii]*sum1[ii];//BODYCOM0 
   sum2[ii] = mhu_a[11];//BODYCOM0 
   sum2[ii] = mhu_a[10]   + VM[ii]*sum2[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[9];//BODYCOM0 
   sum3[ii] = tauRdt_a[8] + VM[ii]*sum3[ii];//BODYCOM0 
   sum4[ii] = tauR_a[16]; //BODYCOM0 
   sum4[ii] = tauR_a[15]  + VM[ii]*sum4[ii]; //BODYCOM0 
*///BODYCOM0 
//BODYCOM0 
  vdt v_mhu_A1;//BODYCOM0 
  vdt v_mhu_A2;//BODYCOM0 
  vdt v_mhu_B1;//BODYCOM0 
  vdt v_mhu_B2;//BODYCOM0 
  vdt v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[8]);//BODYCOM0 
  vdt v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[9]);//BODYCOM0 
  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[15]);//BODYCOM0 
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[16]);//BODYCOM0 
//BODYCOM0 
  vdt v_sum1a = vec_madd(v_xa, v_mhu_a6,  v_mhu_a5);//BODYCOM0 
  vdt v_sum2a = vec_madd(v_xa, v_mhu_a11, v_mhu_a10);//BODYCOM0 
  vdt v_sum3a = vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);//BODYCOM0 
  vdt v_sum4a = vec_madd(v_xa, v_tauR_D2, v_tauR_D1);//BODYCOM0 
//BODYCOM0 
//2//BODYCOM0 
/*//BODYCOM0 
   sum1[ii] = mhu_a[4]    + VM[ii]*sum1[ii];//BODYCOM0 
   sum1[ii] = mhu_a[3]    + VM[ii]*sum1[ii];//BODYCOM0 
   sum2[ii] = mhu_a[9]    + VM[ii]*sum2[ii];//BODYCOM0 
   sum2[ii] = mhu_a[8]    + VM[ii]*sum2[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[7] + VM[ii]*sum3[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[6] + VM[ii]*sum3[ii];//BODYCOM0 
   sum4[ii] = tauR_a[14]  + VM[ii]*sum4[ii]; //BODYCOM0 
   sum4[ii] = tauR_a[13]  + VM[ii]*sum4[ii]; //BODYCOM0 
*///BODYCOM0 
   v_mhu_A1     = VEC_LDS(0, &mhu_a[3]);//BODYCOM0 
   v_mhu_A2     = VEC_LDS(0, &mhu_a[4]);//BODYCOM0 
   v_mhu_B1     = VEC_LDS(0, &mhu_a[8]);//BODYCOM0 
   v_mhu_B2     = VEC_LDS(0, &mhu_a[9]);//BODYCOM0 
   v_tauRdt_C1  = VEC_LDS(0, &tauRdt_a[6]);//BODYCOM0 
   v_tauRdt_C2  = VEC_LDS(0, &tauRdt_a[7]);//BODYCOM0 
   v_tauR_D1    = VEC_LDS(0, &tauR_a[13]);//BODYCOM0 
   v_tauR_D2    = VEC_LDS(0, &tauR_a[14]);//BODYCOM0 
//BODYCOM0 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM0 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM0 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM0 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM0 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM0 
//BODYCOM0 
//3//BODYCOM0 
/*//BODYCOM0 
   sum1[ii] = mhu_a[2]    + VM[ii]*sum1[ii];//BODYCOM0 
   sum1[ii] = mhu_a[1]    + VM[ii]*sum1[ii];//BODYCOM0 
   sum2[ii] = mhu_a[7]    + VM[ii]*sum2[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[5] + VM[ii]*sum3[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[4] + VM[ii]*sum3[ii];//BODYCOM0 
   sum4[ii] = tauR_a[12]  + VM[ii]*sum4[ii]; //BODYCOM0 
   sum4[ii] = tauR_a[11]  + VM[ii]*sum4[ii]; //BODYCOM0 
*///BODYCOM0 
//BODYCOM0 
   v_mhu_A1    = VEC_LDS(0, &mhu_a[1]);//BODYCOM0 
   v_mhu_A2    = VEC_LDS(0, &mhu_a[2]);//BODYCOM0 
   v_mhu_B1    = VEC_LDS(0, &mhu_a[7]);//BODYCOM0 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[4]);//BODYCOM0 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[5]);//BODYCOM0 
   v_tauR_D1   = VEC_LDS(0, &tauR_a[11]);//BODYCOM0 
   v_tauR_D2   = VEC_LDS(0, &tauR_a[12]);//BODYCOM0 
//BODYCOM0 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM0 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM0 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM0 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM0 
//BODYCOM0 
//4//BODYCOM0 
/*//BODYCOM0 
///BODYCOM0 
   sum1[ii] = mhu_a[0]    + VM[ii]*sum1[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[3] + VM[ii]*sum3[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[2] + VM[ii]*sum3[ii];//BODYCOM0 
   sum4[ii] = tauR_a[10]  + VM[ii]*sum4[ii]; //BODYCOM0 
*///BODYCOM0 
   v_mhu_A1    = VEC_LDS(0, &mhu_a[0]);//BODYCOM0 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[2]);//BODYCOM0 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[3]);//BODYCOM0 
   v_tauR_D1   = VEC_LDS(0, &tauR_a[10]);//BODYCOM0 
//BODYCOM0 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM0 
//BODYCOM0 
//5//BODYCOM0 
/*//BODYCOM0 
//BODYCOM0 
   sum3[ii] = tauRdt_a[1]  + VM[ii]*sum3[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[0]  + VM[ii]*sum3[ii];//BODYCOM0 
*///BODYCOM0 
//BODYCOM0 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[0]);//BODYCOM0 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[1]);//BODYCOM0 
//BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
//BODYCOM0 
/*//BODYCOM0 
   double tauRdt= sum3/sum4; //BODYCOM0 
   double mhu= sum1/sum2;//BODYCOM0 
//BODYCOM0 
*///BODYCOM0 
//BODYCOM0 
/*//BODYCOM0 
   tauRdt[ii] = sum3[ii]/sum4[ii]; //BODYCOM0 
   mhu[ii]= sum1[ii]/sum2[ii];//BODYCOM0 
//BODYCOM0 
   g[ii] +=  mhu[ii]*tauRdt[ii] - g[ii]*tauRdt[ii];//BODYCOM0 
*///BODYCOM0 
//BODYCOM0 
// Still need to reuse the register I have above ... the code will look like gibberish because//BODYCOM0 
// the targets won't make sense.//BODYCOM0 
// A//BODYCOM0 
   v_xa      = vec_swdiv_nochk(v_sum3a, v_sum4a); //BODYCOM0 
   v_mhu_A2  = vec_swdiv_nochk(v_sum1a, v_sum2a); //BODYCOM0 
   v_mhu_B2  = vec_ld(0, &g[ii]);  //BODYCOM0 
   v_tauR_D2 = vec_sub(v_mhu_A2, v_mhu_B2);  //BODYCOM0 
   v_mhu_B2  = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);//BODYCOM0 

   
   //vec_st(v_mhu_B2, 0, &g[ii]);
   switch (misalignment) {
   case 1:
     g[ii+0] = v_mhu_B2 get [1];
     g[ii+1] = v_mhu_B2 get [2];
     g[ii+2] = v_mhu_B2 get [3];
     break;

   case 2:
     g[ii+0] = v_mhu_B2 get [2];
     g[ii+1] = v_mhu_B2 get [3];
     break;

   case 3:
     g[ii+0] = v_mhu_B2 get [3];
     break;

     default:
    //println("BAD MISALIGNMENT: %d\n", misalignment); fflush(stdout);
     assert(0);
   }

   ii = 4 - misalignment;
 }

 for (ii;ii<TTROUND(nCells,8);ii+=8)
 { 
   vdt v_xa   = vec_ld(0, &VM[ii]);  
   vdt v_xb   = vec_ld(0, &VM[ii+4]);  
// 1
/*
   sum1[ii] = mhu_a[6];
   sum1[ii] = mhu_a[5]    + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[11];
   sum2[ii] = mhu_a[10]   + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[9];
   sum3[ii] = tauRdt_a[8] + VM[ii]*sum3[ii];
   sum4[ii] = tauR_a[16]; 
   sum4[ii] = tauR_a[15]  + VM[ii]*sum4[ii]; 
*/

  vdt v_mhu_A1;
  vdt v_mhu_A2;
  vdt v_mhu_B1;
  vdt v_mhu_B2;
  vdt v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[8]);
  vdt v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[9]);
  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[15]);
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[16]);

  vdt v_sum1a = vec_madd(v_xa, v_mhu_a6,  v_mhu_a5);
  vdt v_sum1b = vec_madd(v_xb, v_mhu_a6,  v_mhu_a5);
  vdt v_sum2a = vec_madd(v_xa, v_mhu_a11, v_mhu_a10);
  vdt v_sum2b = vec_madd(v_xb, v_mhu_a11, v_mhu_a10);
  vdt v_sum3a = vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);
  vdt v_sum3b = vec_madd(v_xb, v_tauRdt_C2, v_tauRdt_C1);
  vdt v_sum4a = vec_madd(v_xa, v_tauR_D2, v_tauR_D1);
  vdt v_sum4b = vec_madd(v_xb, v_tauR_D2, v_tauR_D1);

//2
/*
   sum1[ii] = mhu_a[4]    + VM[ii]*sum1[ii];
   sum1[ii] = mhu_a[3]    + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[9]    + VM[ii]*sum2[ii];
   sum2[ii] = mhu_a[8]    + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[7] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[6] + VM[ii]*sum3[ii];
   sum4[ii] = tauR_a[14]  + VM[ii]*sum4[ii]; 
   sum4[ii] = tauR_a[13]  + VM[ii]*sum4[ii]; 
*/
   v_mhu_A1     = VEC_LDS(0, &mhu_a[3]);
   v_mhu_A2     = VEC_LDS(0, &mhu_a[4]);
   v_mhu_B1     = VEC_LDS(0, &mhu_a[8]);
   v_mhu_B2     = VEC_LDS(0, &mhu_a[9]);
   v_tauRdt_C1  = VEC_LDS(0, &tauRdt_a[6]);
   v_tauRdt_C2  = VEC_LDS(0, &tauRdt_a[7]);
   v_tauR_D1    = VEC_LDS(0, &tauR_a[13]);
   v_tauR_D2    = VEC_LDS(0, &tauR_a[14]);

   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);

   v_sum1b = vec_madd(v_xb, v_sum1b, v_mhu_A2);
   v_sum1b = vec_madd(v_xb, v_sum1b, v_mhu_A1);
   v_sum2b = vec_madd(v_xb, v_sum2b, v_mhu_B2);
   v_sum2b = vec_madd(v_xb, v_sum2b, v_mhu_B1); 
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   v_sum4b = vec_madd(v_xb, v_sum4b, v_tauR_D2);
   v_sum4b = vec_madd(v_xb, v_sum4b, v_tauR_D1);

//3
/*
   sum1[ii] = mhu_a[2]    + VM[ii]*sum1[ii];
   sum1[ii] = mhu_a[1]    + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[7]    + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[5] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[4] + VM[ii]*sum3[ii];
   sum4[ii] = tauR_a[12]  + VM[ii]*sum4[ii]; 
   sum4[ii] = tauR_a[11]  + VM[ii]*sum4[ii]; 
*/

   v_mhu_A1    = VEC_LDS(0, &mhu_a[1]);
   v_mhu_A2    = VEC_LDS(0, &mhu_a[2]);
   v_mhu_B1    = VEC_LDS(0, &mhu_a[7]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[4]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[5]);
   v_tauR_D1   = VEC_LDS(0, &tauR_a[11]);
   v_tauR_D2   = VEC_LDS(0, &tauR_a[12]);

   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);

   v_sum1b = vec_madd(v_xb, v_sum1b, v_mhu_A2);
   v_sum1b = vec_madd(v_xb, v_sum1b, v_mhu_A1);
   v_sum2b = vec_madd(v_xb, v_sum2b, v_mhu_B1); 
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   v_sum4b = vec_madd(v_xb, v_sum4b, v_tauR_D2);
   v_sum4b = vec_madd(v_xb, v_sum4b, v_tauR_D1);

//4
/*
/
   sum1[ii] = mhu_a[0]    + VM[ii]*sum1[ii];
   sum3[ii] = tauRdt_a[3] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[2] + VM[ii]*sum3[ii];
   sum4[ii] = tauR_a[10]  + VM[ii]*sum4[ii]; 
*/
   v_mhu_A1    = VEC_LDS(0, &mhu_a[0]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[2]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[3]);
   v_tauR_D1   = VEC_LDS(0, &tauR_a[10]);

   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);

   v_sum1b = vec_madd(v_xb, v_sum1b, v_mhu_A1);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   v_sum4b = vec_madd(v_xb, v_sum4b, v_tauR_D1);

//5
/*

   sum3[ii] = tauRdt_a[1]  + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[0]  + VM[ii]*sum3[ii];
*/

   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[0]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[1]);

   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);

   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
 

/*
   double tauRdt= sum3/sum4; 
   double mhu= sum1/sum2;

*/

/*
   tauRdt[ii] = sum3[ii]/sum4[ii]; 
   mhu[ii]= sum1[ii]/sum2[ii];

   g[ii] +=  mhu[ii]*tauRdt[ii] - g[ii]*tauRdt[ii];
*/

// Still need to reuse the register I have above ... the code will look like gibberish because
// the targets won't make sense.
// A
   v_xa      = vec_swdiv_nochk(v_sum3a, v_sum4a); 
   v_mhu_A2  = vec_swdiv_nochk(v_sum1a, v_sum2a); 
   v_mhu_B2  = vec_ld(0, &g[ii]);  
   v_tauR_D2 = vec_sub(v_mhu_A2, v_mhu_B2);  
   v_mhu_B2  = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);
   vec_st(v_mhu_B2, 0, &g[ii]);

// B
   v_xb      = vec_swdiv_nochk(v_sum3b, v_sum4b); 
   v_mhu_A1  = vec_swdiv_nochk(v_sum1b, v_sum2b); 
   v_mhu_B1  = vec_ld(0, &g[ii+4]);  
   v_tauR_D1 = vec_sub(v_mhu_A1, v_mhu_B1);  
   v_mhu_B1  = vec_madd(v_tauR_D1, v_xb, v_mhu_B1);
   vec_st(v_mhu_B1, 0, &g[ii+4]);



  } 
 for (ii;ii<TTROUND(nCells,4);ii+=4)
 { 
 //BODYSTART
   vdt v_xa   = vec_ld(0, &VM[ii]);  

// 1
/*
   sum1[ii] = mhu_a[6];
   sum1[ii] = mhu_a[5]    + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[11];
   sum2[ii] = mhu_a[10]   + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[9];
   sum3[ii] = tauRdt_a[8] + VM[ii]*sum3[ii];
   sum4[ii] = tauR_a[16]; 
   sum4[ii] = tauR_a[15]  + VM[ii]*sum4[ii]; 
*/

  vdt v_mhu_A1;
  vdt v_mhu_A2;
  vdt v_mhu_B1;
  vdt v_mhu_B2;
  vdt v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[8]);
  vdt v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[9]);
  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[15]);
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[16]);

  vdt v_sum1a = vec_madd(v_xa, v_mhu_a6,  v_mhu_a5);
  vdt v_sum2a = vec_madd(v_xa, v_mhu_a11, v_mhu_a10);
  vdt v_sum3a = vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);
  vdt v_sum4a = vec_madd(v_xa, v_tauR_D2, v_tauR_D1);

//2
/*
   sum1[ii] = mhu_a[4]    + VM[ii]*sum1[ii];
   sum1[ii] = mhu_a[3]    + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[9]    + VM[ii]*sum2[ii];
   sum2[ii] = mhu_a[8]    + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[7] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[6] + VM[ii]*sum3[ii];
   sum4[ii] = tauR_a[14]  + VM[ii]*sum4[ii]; 
   sum4[ii] = tauR_a[13]  + VM[ii]*sum4[ii]; 
*/
   v_mhu_A1     = VEC_LDS(0, &mhu_a[3]);
   v_mhu_A2     = VEC_LDS(0, &mhu_a[4]);
   v_mhu_B1     = VEC_LDS(0, &mhu_a[8]);
   v_mhu_B2     = VEC_LDS(0, &mhu_a[9]);
   v_tauRdt_C1  = VEC_LDS(0, &tauRdt_a[6]);
   v_tauRdt_C2  = VEC_LDS(0, &tauRdt_a[7]);
   v_tauR_D1    = VEC_LDS(0, &tauR_a[13]);
   v_tauR_D2    = VEC_LDS(0, &tauR_a[14]);

   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);

//3
/*
   sum1[ii] = mhu_a[2]    + VM[ii]*sum1[ii];
   sum1[ii] = mhu_a[1]    + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[7]    + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[5] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[4] + VM[ii]*sum3[ii];
   sum4[ii] = tauR_a[12]  + VM[ii]*sum4[ii]; 
   sum4[ii] = tauR_a[11]  + VM[ii]*sum4[ii]; 
*/

   v_mhu_A1    = VEC_LDS(0, &mhu_a[1]);
   v_mhu_A2    = VEC_LDS(0, &mhu_a[2]);
   v_mhu_B1    = VEC_LDS(0, &mhu_a[7]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[4]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[5]);
   v_tauR_D1   = VEC_LDS(0, &tauR_a[11]);
   v_tauR_D2   = VEC_LDS(0, &tauR_a[12]);

   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);

//4
/*
/
   sum1[ii] = mhu_a[0]    + VM[ii]*sum1[ii];
   sum3[ii] = tauRdt_a[3] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[2] + VM[ii]*sum3[ii];
   sum4[ii] = tauR_a[10]  + VM[ii]*sum4[ii]; 
*/
   v_mhu_A1    = VEC_LDS(0, &mhu_a[0]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[2]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[3]);
   v_tauR_D1   = VEC_LDS(0, &tauR_a[10]);

   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);

//5
/*

   sum3[ii] = tauRdt_a[1]  + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[0]  + VM[ii]*sum3[ii];
*/

   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[0]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[1]);

   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);

/*
   double tauRdt= sum3/sum4; 
   double mhu= sum1/sum2;

*/

/*
   tauRdt[ii] = sum3[ii]/sum4[ii]; 
   mhu[ii]= sum1[ii]/sum2[ii];

   g[ii] +=  mhu[ii]*tauRdt[ii] - g[ii]*tauRdt[ii];
*/

// Still need to reuse the register I have above ... the code will look like gibberish because
// the targets won't make sense.
// A
   v_xa      = vec_swdiv_nochk(v_sum3a, v_sum4a); 
   v_mhu_A2  = vec_swdiv_nochk(v_sum1a, v_sum2a); 
   v_mhu_B2  = vec_ld(0, &g[ii]);  
   v_tauR_D2 = vec_sub(v_mhu_A2, v_mhu_B2);  
   v_mhu_B2  = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);
   vec_st(v_mhu_B2, 0, &g[ii]); //BODYEND
  } 

 // the elements in the end of the vector should be ignored
 int leftover = nCells - ii;
 if (leftover) { 
      vdt v_xa   = vec_ld(0, &VM[ii]);  //BODYCOM1 
//BODYCOM1 
// 1//BODYCOM1 
/*//BODYCOM1 
   sum1[ii] = mhu_a[6];//BODYCOM1 
   sum1[ii] = mhu_a[5]    + VM[ii]*sum1[ii];//BODYCOM1 
   sum2[ii] = mhu_a[11];//BODYCOM1 
   sum2[ii] = mhu_a[10]   + VM[ii]*sum2[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[9];//BODYCOM1 
   sum3[ii] = tauRdt_a[8] + VM[ii]*sum3[ii];//BODYCOM1 
   sum4[ii] = tauR_a[16]; //BODYCOM1 
   sum4[ii] = tauR_a[15]  + VM[ii]*sum4[ii]; //BODYCOM1 
*///BODYCOM1 
//BODYCOM1 
  vdt v_mhu_A1;//BODYCOM1 
  vdt v_mhu_A2;//BODYCOM1 
  vdt v_mhu_B1;//BODYCOM1 
  vdt v_mhu_B2;//BODYCOM1 
  vdt v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[8]);//BODYCOM1 
  vdt v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[9]);//BODYCOM1 
  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[15]);//BODYCOM1 
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[16]);//BODYCOM1 
//BODYCOM1 
  vdt v_sum1a = vec_madd(v_xa, v_mhu_a6,  v_mhu_a5);//BODYCOM1 
  vdt v_sum2a = vec_madd(v_xa, v_mhu_a11, v_mhu_a10);//BODYCOM1 
  vdt v_sum3a = vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);//BODYCOM1 
  vdt v_sum4a = vec_madd(v_xa, v_tauR_D2, v_tauR_D1);//BODYCOM1 
//BODYCOM1 
//2//BODYCOM1 
/*//BODYCOM1 
   sum1[ii] = mhu_a[4]    + VM[ii]*sum1[ii];//BODYCOM1 
   sum1[ii] = mhu_a[3]    + VM[ii]*sum1[ii];//BODYCOM1 
   sum2[ii] = mhu_a[9]    + VM[ii]*sum2[ii];//BODYCOM1 
   sum2[ii] = mhu_a[8]    + VM[ii]*sum2[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[7] + VM[ii]*sum3[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[6] + VM[ii]*sum3[ii];//BODYCOM1 
   sum4[ii] = tauR_a[14]  + VM[ii]*sum4[ii]; //BODYCOM1 
   sum4[ii] = tauR_a[13]  + VM[ii]*sum4[ii]; //BODYCOM1 
*///BODYCOM1 
   v_mhu_A1     = VEC_LDS(0, &mhu_a[3]);//BODYCOM1 
   v_mhu_A2     = VEC_LDS(0, &mhu_a[4]);//BODYCOM1 
   v_mhu_B1     = VEC_LDS(0, &mhu_a[8]);//BODYCOM1 
   v_mhu_B2     = VEC_LDS(0, &mhu_a[9]);//BODYCOM1 
   v_tauRdt_C1  = VEC_LDS(0, &tauRdt_a[6]);//BODYCOM1 
   v_tauRdt_C2  = VEC_LDS(0, &tauRdt_a[7]);//BODYCOM1 
   v_tauR_D1    = VEC_LDS(0, &tauR_a[13]);//BODYCOM1 
   v_tauR_D2    = VEC_LDS(0, &tauR_a[14]);//BODYCOM1 
//BODYCOM1 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM1 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM1 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM1 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM1 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM1 
//BODYCOM1 
//3//BODYCOM1 
/*//BODYCOM1 
   sum1[ii] = mhu_a[2]    + VM[ii]*sum1[ii];//BODYCOM1 
   sum1[ii] = mhu_a[1]    + VM[ii]*sum1[ii];//BODYCOM1 
   sum2[ii] = mhu_a[7]    + VM[ii]*sum2[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[5] + VM[ii]*sum3[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[4] + VM[ii]*sum3[ii];//BODYCOM1 
   sum4[ii] = tauR_a[12]  + VM[ii]*sum4[ii]; //BODYCOM1 
   sum4[ii] = tauR_a[11]  + VM[ii]*sum4[ii]; //BODYCOM1 
*///BODYCOM1 
//BODYCOM1 
   v_mhu_A1    = VEC_LDS(0, &mhu_a[1]);//BODYCOM1 
   v_mhu_A2    = VEC_LDS(0, &mhu_a[2]);//BODYCOM1 
   v_mhu_B1    = VEC_LDS(0, &mhu_a[7]);//BODYCOM1 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[4]);//BODYCOM1 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[5]);//BODYCOM1 
   v_tauR_D1   = VEC_LDS(0, &tauR_a[11]);//BODYCOM1 
   v_tauR_D2   = VEC_LDS(0, &tauR_a[12]);//BODYCOM1 
//BODYCOM1 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM1 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM1 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM1 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM1 
//BODYCOM1 
//4//BODYCOM1 
/*//BODYCOM1 
///BODYCOM1 
   sum1[ii] = mhu_a[0]    + VM[ii]*sum1[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[3] + VM[ii]*sum3[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[2] + VM[ii]*sum3[ii];//BODYCOM1 
   sum4[ii] = tauR_a[10]  + VM[ii]*sum4[ii]; //BODYCOM1 
*///BODYCOM1 
   v_mhu_A1    = VEC_LDS(0, &mhu_a[0]);//BODYCOM1 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[2]);//BODYCOM1 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[3]);//BODYCOM1 
   v_tauR_D1   = VEC_LDS(0, &tauR_a[10]);//BODYCOM1 
//BODYCOM1 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM1 
//BODYCOM1 
//5//BODYCOM1 
/*//BODYCOM1 
//BODYCOM1 
   sum3[ii] = tauRdt_a[1]  + VM[ii]*sum3[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[0]  + VM[ii]*sum3[ii];//BODYCOM1 
*///BODYCOM1 
//BODYCOM1 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[0]);//BODYCOM1 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[1]);//BODYCOM1 
//BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
//BODYCOM1 
/*//BODYCOM1 
   double tauRdt= sum3/sum4; //BODYCOM1 
   double mhu= sum1/sum2;//BODYCOM1 
//BODYCOM1 
*///BODYCOM1 
//BODYCOM1 
/*//BODYCOM1 
   tauRdt[ii] = sum3[ii]/sum4[ii]; //BODYCOM1 
   mhu[ii]= sum1[ii]/sum2[ii];//BODYCOM1 
//BODYCOM1 
   g[ii] +=  mhu[ii]*tauRdt[ii] - g[ii]*tauRdt[ii];//BODYCOM1 
*///BODYCOM1 
//BODYCOM1 
// Still need to reuse the register I have above ... the code will look like gibberish because//BODYCOM1 
// the targets won't make sense.//BODYCOM1 
// A//BODYCOM1 
   v_xa      = vec_swdiv_nochk(v_sum3a, v_sum4a); //BODYCOM1 
   v_mhu_A2  = vec_swdiv_nochk(v_sum1a, v_sum2a); //BODYCOM1 
   v_mhu_B2  = vec_ld(0, &g[ii]);  //BODYCOM1 
   v_tauR_D2 = vec_sub(v_mhu_A2, v_mhu_B2);  //BODYCOM1 
   v_mhu_B2  = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);//BODYCOM1 

   
   //vec_st(v_mhu_B2, 0, &g[ii]);
   switch (leftover) {
   case 1:
     g[ii+0] = v_mhu_B2 get [0];
     break;

   case 2:
     g[ii+0] = v_mhu_B2 get [0];
     g[ii+1] = v_mhu_B2 get [1];
     break;

   case 3:
     g[ii+0] = v_mhu_B2 get [0];
     g[ii+1] = v_mhu_B2 get [1];
     g[ii+2] = v_mhu_B2 get [2];
     break;

     default:
    //println("BAD LEFTOVER: %d\n", leftover); fflush(stdout);
    assert(0);
   }
 }

}





void update_fGate_v1(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a)
{
  //update_fGate(dt, nCells, VM, g, mhu_a, tauR_a); return;
  typedef vector4double vdt;

int gateIndex=8;
int mhu_l=5 ;
int mhu_m= 8;
int tauR_l=13;
int tauR_m=11;
   
 double  tauRdt_a[11];


         tauRdt_a[0]  = tauR_a[0]*dt;
         tauRdt_a[1]  = tauR_a[1]*dt;
         tauRdt_a[2]  = tauR_a[2]*dt;
         tauRdt_a[3]  = tauR_a[3]*dt;
         tauRdt_a[4]  = tauR_a[4]*dt;
         tauRdt_a[5]  = tauR_a[5]*dt;
         tauRdt_a[6]  = tauR_a[6]*dt; 
         tauRdt_a[7]  = tauR_a[7]*dt;
         tauRdt_a[8]  = tauR_a[8]*dt;
         tauRdt_a[9]  = tauR_a[9]*dt;
         tauRdt_a[10] = tauR_a[10]*dt;

 int ii=0;

 // 0 aligned
 // 1: load produces (JUNK, &VM[0], &VM[1], &VM[2])
 // 2: load produces (JUNK, JUNK, &VM[0], &VM[1])
 // 3: load produces (JUNK, JUNK, JUNK, &VM[0])
 int misalignment = ((long int)&VM[0] % 32) / 8;
 assert(misalignment == 0);

 if (misalignment) { 
      vdt v_xa   = vec_ld(0, &VM[ii]);  //BODYCOM0 
//BODYCOM0 
// 1//BODYCOM0 
  vdt v_mhu_A1  = VEC_LDS(0, &mhu_a[6]);//BODYCOM0 
  vdt v_mhu_A2  = VEC_LDS(0, &mhu_a[7]);//BODYCOM0 
  vdt v_mhu_B1 = VEC_LDS(0, &mhu_a[11]);//BODYCOM0 
  vdt v_mhu_B2 = VEC_LDS(0, &mhu_a[12]);//BODYCOM0 
  vdt v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[9]);//BODYCOM0 
  vdt v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[10]);//BODYCOM0 
  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[22]);//BODYCOM0 
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[23]);//BODYCOM0 
//BODYCOM0 
  vdt v_sum1a =  vec_madd(v_xa, v_mhu_A2, v_mhu_A1);//BODYCOM0 
   vdt v_sum2a =  vec_madd(v_xa, v_mhu_B2, v_mhu_B1);//BODYCOM0 
   vdt v_sum3a =  vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);//BODYCOM0 
   vdt v_sum4a =  vec_madd(v_xa, v_tauR_D2, v_tauR_D1);//BODYCOM0 
// 2//BODYCOM0 
   v_mhu_A1  = VEC_LDS(0, &mhu_a[4]);//BODYCOM0 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[5]);//BODYCOM0 
   v_mhu_B1 = VEC_LDS(0, &mhu_a[9]);//BODYCOM0 
   v_mhu_B2 = VEC_LDS(0, &mhu_a[10]);//BODYCOM0 
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[7]);//BODYCOM0 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[8]);//BODYCOM0 
  v_tauR_D1   = VEC_LDS(0, &tauR_a[20]);//BODYCOM0 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[21]);//BODYCOM0 
//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM0 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM0 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM0 
// 3//BODYCOM0 
//BODYCOM0 
   v_mhu_A1  = VEC_LDS(0, &mhu_a[2]);//BODYCOM0 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[3]);//BODYCOM0 
   v_mhu_B2 = VEC_LDS(0, &mhu_a[8]);//BODYCOM0 
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[5]);//BODYCOM0 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[6]);//BODYCOM0 
  v_tauR_D1   = VEC_LDS(0, &tauR_a[18]);//BODYCOM0 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[19]);//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM0 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM0 
// 4//BODYCOM0 
   v_mhu_A1  = VEC_LDS(0, &mhu_a[0]);//BODYCOM0 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[1]);//BODYCOM0 
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[3]);//BODYCOM0 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[4]);//BODYCOM0 
  v_tauR_D1   = VEC_LDS(0, &tauR_a[16]);//BODYCOM0 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[17]);//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM0 
// 5//BODYCOM0 
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[1]);//BODYCOM0 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[2]);//BODYCOM0 
  v_tauR_D1   = VEC_LDS(0, &tauR_a[14]);//BODYCOM0 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[15]);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM0 
// 6//BODYCOM0 
   v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[0]);//BODYCOM0 
   v_tauR_D1   = VEC_LDS(0, &tauR_a[12]);//BODYCOM0 
   v_tauR_D2   = VEC_LDS(0, &tauR_a[13]);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2); // Was missing//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM0 
// 7//BODYCOM0 
   v_tauR_D2   = VEC_LDS(0, &tauR_a[11]);//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM0 
//BODYCOM0 
// A//BODYCOM0 
   vdt temp1a   = vec_swdiv_nochk(v_sum3a,v_sum4a); //BODYCOM0 
   vdt temp2a   = vec_swdiv_nochk(v_sum1a,v_sum2a); //BODYCOM0 
   vdt temp3a   = vec_ld(0, &g[ii]);  //BODYCOM0 
   temp2a   = vec_sub(temp2a, temp3a);  //BODYCOM0 
   temp3a   = vec_madd(temp2a, temp1a, temp3a);//BODYCOM0 

   
   //vec_st(temp3a, 0, &g[ii]);
   switch (misalignment) {
   case 1:
     g[ii+0] = temp3a get [1];
     g[ii+1] = temp3a get [2];
     g[ii+2] = temp3a get [3];
     break;

   case 2:
     g[ii+0] = temp3a get [2];
     g[ii+1] = temp3a get [3];
     break;

   case 3:
     g[ii+0] = temp3a get [3];
     break;

     default:
    //println("BAD MISALIGNMENT: %d\n", misalignment); fflush(stdout);
     assert(0);
   }

   ii = 4 - misalignment;
 }

 //for (ii=0;ii<0;ii+=16)
 for (ii;ii<TTROUND(nCells,12);ii+=12)
 { 
   vdt v_xa   = vec_ld(0, &VM[ii]);  
   vdt v_xb   = vec_ld(0, &VM[ii+4]);  
   vdt v_xc   = vec_ld(0, &VM[ii+8]);  
// 1
  vdt v_mhu_A1  = VEC_LDS(0, &mhu_a[6]);
  vdt v_mhu_A2  = VEC_LDS(0, &mhu_a[7]);
  vdt v_mhu_B1 = VEC_LDS(0, &mhu_a[11]);
  vdt v_mhu_B2 = VEC_LDS(0, &mhu_a[12]);
  vdt v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[9]);
  vdt v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[10]);
  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[22]);
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[23]);

   vdt v_sum1a =  vec_madd(v_xa, v_mhu_A2, v_mhu_A1);
   vdt v_sum1b =  vec_madd(v_xb, v_mhu_A2, v_mhu_A1);
   vdt v_sum1c =  vec_madd(v_xc, v_mhu_A2, v_mhu_A1);

   vdt v_sum2a =  vec_madd(v_xa, v_mhu_B2, v_mhu_B1);
   vdt v_sum2b =  vec_madd(v_xb, v_mhu_B2, v_mhu_B1);
   vdt v_sum2c =  vec_madd(v_xc, v_mhu_B2, v_mhu_B1);

   vdt v_sum3a =  vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);
   vdt v_sum3b =  vec_madd(v_xb, v_tauRdt_C2, v_tauRdt_C1);
   vdt v_sum3c =  vec_madd(v_xc, v_tauRdt_C2, v_tauRdt_C1);

   vdt v_sum4a =  vec_madd(v_xa, v_tauR_D2, v_tauR_D1);
   vdt v_sum4b =  vec_madd(v_xb, v_tauR_D2, v_tauR_D1);
   vdt v_sum4c =  vec_madd(v_xc, v_tauR_D2, v_tauR_D1);
// 2
   v_mhu_A1  = VEC_LDS(0, &mhu_a[4]);
   v_mhu_A2  = VEC_LDS(0, &mhu_a[5]);
   v_mhu_B1 = VEC_LDS(0, &mhu_a[9]);
   v_mhu_B2 = VEC_LDS(0, &mhu_a[10]);
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[7]);
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[8]);
  v_tauR_D1   = VEC_LDS(0, &tauR_a[20]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[21]);

   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A2);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A1);
   v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B2);
   v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B1); 
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D2);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D1);
   v_sum1c =  vec_madd(v_xc, v_sum1c, v_mhu_A2);
   v_sum1c =  vec_madd(v_xc, v_sum1c, v_mhu_A1);
   v_sum2c =  vec_madd(v_xc, v_sum2c, v_mhu_B2);
   v_sum2c =  vec_madd(v_xc, v_sum2c, v_mhu_B1); 
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C2);
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C1);
   v_sum4c =  vec_madd(v_xc, v_sum4c, v_tauR_D2);
   v_sum4c =  vec_madd(v_xc, v_sum4c, v_tauR_D1);
// 3

   v_mhu_A1  = VEC_LDS(0, &mhu_a[2]);
   v_mhu_A2  = VEC_LDS(0, &mhu_a[3]);
   v_mhu_B2 = VEC_LDS(0, &mhu_a[8]);
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[5]);
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[6]);
  v_tauR_D1   = VEC_LDS(0, &tauR_a[18]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[19]);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A2);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A1);
   v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B2);
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D2);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D1);
   v_sum1c =  vec_madd(v_xc, v_sum1c, v_mhu_A2);
   v_sum1c =  vec_madd(v_xc, v_sum1c, v_mhu_A1);
   v_sum2c =  vec_madd(v_xc, v_sum2c, v_mhu_B2);
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C2);
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C1);
   v_sum4c =  vec_madd(v_xc, v_sum4c, v_tauR_D2);
   v_sum4c =  vec_madd(v_xc, v_sum4c, v_tauR_D1);
// 4

   v_mhu_A1  = VEC_LDS(0, &mhu_a[0]);
   v_mhu_A2  = VEC_LDS(0, &mhu_a[1]);
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[3]);
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[4]);
  v_tauR_D1   = VEC_LDS(0, &tauR_a[16]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[17]);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A2);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A1);
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D2);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D1);
   v_sum1c =  vec_madd(v_xc, v_sum1c, v_mhu_A2);
   v_sum1c =  vec_madd(v_xc, v_sum1c, v_mhu_A1);
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C2);
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C1);
   v_sum4c =  vec_madd(v_xc, v_sum4c, v_tauR_D2);
   v_sum4c =  vec_madd(v_xc, v_sum4c, v_tauR_D1);
// 5
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[1]);
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[2]);
  v_tauR_D1   = VEC_LDS(0, &tauR_a[14]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[15]);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D2);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D1);
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C2);
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C1);
   v_sum4c =  vec_madd(v_xc, v_sum4c, v_tauR_D2);
   v_sum4c =  vec_madd(v_xc, v_sum4c, v_tauR_D1);
// 6
   v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[0]);
   v_tauR_D1   = VEC_LDS(0, &tauR_a[12]);
   v_tauR_D2   = VEC_LDS(0, &tauR_a[13]);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C2);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D2);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D1);
   v_sum4c =  vec_madd(v_xc, v_sum4c, v_tauR_D2);
   v_sum4c =  vec_madd(v_xc, v_sum4c, v_tauR_D1);
// 7
   v_tauR_D2   = VEC_LDS(0, &tauR_a[11]);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D2);
   v_sum4c =  vec_madd(v_xc, v_sum4c, v_tauR_D2);

// A
   vdt temp1a   = vec_swdiv_nochk(v_sum3a,v_sum4a); 
   vdt temp2a   = vec_swdiv_nochk(v_sum1a,v_sum2a); 
   vdt temp3a   = vec_ld(0, &g[ii]);  
   temp2a   = vec_sub(temp2a, temp3a);  
   temp3a   = vec_madd(temp2a, temp1a, temp3a);
   vec_st(temp3a, 0, &g[ii]);
   /* This can also be implemented as the following, which has a higher flop rate and is more straightforward
   temp1a  = vec_swdiv_nochk(v_sum3a,v_sum4a); 
   temp2a  = vec_swdiv_nochk(v_sum1a,v_sum2a); 
   temp3a  = vec_ld(0, &g[ii]);
   temp2a   = vec_madd(temp1a, temp2a, temp3a);
   temp3a  = vec_nmadd(temp1a, temp3a, temp2a );
   vec_st(temp3a, 0, &g[ii]); */
// B
   vdt temp1b   = vec_swdiv_nochk(v_sum3b,v_sum4b); 
   vdt temp2b   = vec_swdiv_nochk(v_sum1b,v_sum2b); 
   vdt temp3b   = vec_ld(0, &g[ii+4]);  
   temp2b   = vec_sub(temp2b, temp3b);  
   temp3b   = vec_madd(temp2b, temp1b, temp3b);
   vec_st(temp3b, 0, &g[ii+4]);
// C
   vdt temp1c   = vec_swdiv_nochk(v_sum3c,v_sum4c); 
   vdt temp2c   = vec_swdiv_nochk(v_sum1c,v_sum2c); 
   vdt temp3c   = vec_ld(0, &g[ii+8]);  
   temp2c   = vec_sub(temp2c, temp3c);  
   temp3c   = vec_madd(temp2c, temp1c, temp3c);
   vec_st(temp3c, 0, &g[ii+8]);
 }
 for (ii;ii<TTROUND(nCells,4);ii+=4)
 { 
 //BODYSTART
   vdt v_xa   = vec_ld(0, &VM[ii]);  

// 1
  vdt v_mhu_A1  = VEC_LDS(0, &mhu_a[6]);
  vdt v_mhu_A2  = VEC_LDS(0, &mhu_a[7]);
  vdt v_mhu_B1 = VEC_LDS(0, &mhu_a[11]);
  vdt v_mhu_B2 = VEC_LDS(0, &mhu_a[12]);
  vdt v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[9]);
  vdt v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[10]);
  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[22]);
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[23]);

  vdt v_sum1a =  vec_madd(v_xa, v_mhu_A2, v_mhu_A1);
   vdt v_sum2a =  vec_madd(v_xa, v_mhu_B2, v_mhu_B1);
   vdt v_sum3a =  vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);
   vdt v_sum4a =  vec_madd(v_xa, v_tauR_D2, v_tauR_D1);
// 2
   v_mhu_A1  = VEC_LDS(0, &mhu_a[4]);
   v_mhu_A2  = VEC_LDS(0, &mhu_a[5]);
   v_mhu_B1 = VEC_LDS(0, &mhu_a[9]);
   v_mhu_B2 = VEC_LDS(0, &mhu_a[10]);
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[7]);
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[8]);
  v_tauR_D1   = VEC_LDS(0, &tauR_a[20]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[21]);

   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);
// 3

   v_mhu_A1  = VEC_LDS(0, &mhu_a[2]);
   v_mhu_A2  = VEC_LDS(0, &mhu_a[3]);
   v_mhu_B2 = VEC_LDS(0, &mhu_a[8]);
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[5]);
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[6]);
  v_tauR_D1   = VEC_LDS(0, &tauR_a[18]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[19]);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);
// 4
   v_mhu_A1  = VEC_LDS(0, &mhu_a[0]);
   v_mhu_A2  = VEC_LDS(0, &mhu_a[1]);
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[3]);
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[4]);
  v_tauR_D1   = VEC_LDS(0, &tauR_a[16]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[17]);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);
// 5
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[1]);
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[2]);
  v_tauR_D1   = VEC_LDS(0, &tauR_a[14]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[15]);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);
// 6
   v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[0]);
   v_tauR_D1   = VEC_LDS(0, &tauR_a[12]);
   v_tauR_D2   = VEC_LDS(0, &tauR_a[13]);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2); // Was missing
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);
// 7
   v_tauR_D2   = VEC_LDS(0, &tauR_a[11]);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);

// A
   vdt temp1a   = vec_swdiv_nochk(v_sum3a,v_sum4a); 
   vdt temp2a   = vec_swdiv_nochk(v_sum1a,v_sum2a); 
   vdt temp3a   = vec_ld(0, &g[ii]);  
   temp2a   = vec_sub(temp2a, temp3a);  
   temp3a   = vec_madd(temp2a, temp1a, temp3a);
   vec_st(temp3a, 0, &g[ii]); //BODYEND
   /* This can also be implemented as the following, which has a higher flop rate and is more straightforward
   temp1a  = vec_swdiv_nochk(v_sum3a,v_sum4a); 
   temp2a  = vec_swdiv_nochk(v_sum1a,v_sum2a); 
   temp3a  = vec_ld(0, &g[ii]);
   temp2a   = vec_madd(temp1a, temp2a, temp3a);
   temp3a  = vec_nmadd(temp1a, temp3a, temp2a );
   vec_st(temp3a, 0, &g[ii]); */
 }

 // the elements in the end of the vector should be ignored
 int leftover = nCells - ii;
 if (leftover) { 
      vdt v_xa   = vec_ld(0, &VM[ii]);  //BODYCOM1 
//BODYCOM1 
// 1//BODYCOM1 
  vdt v_mhu_A1  = VEC_LDS(0, &mhu_a[6]);//BODYCOM1 
  vdt v_mhu_A2  = VEC_LDS(0, &mhu_a[7]);//BODYCOM1 
  vdt v_mhu_B1 = VEC_LDS(0, &mhu_a[11]);//BODYCOM1 
  vdt v_mhu_B2 = VEC_LDS(0, &mhu_a[12]);//BODYCOM1 
  vdt v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[9]);//BODYCOM1 
  vdt v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[10]);//BODYCOM1 
  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[22]);//BODYCOM1 
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[23]);//BODYCOM1 
//BODYCOM1 
  vdt v_sum1a =  vec_madd(v_xa, v_mhu_A2, v_mhu_A1);//BODYCOM1 
   vdt v_sum2a =  vec_madd(v_xa, v_mhu_B2, v_mhu_B1);//BODYCOM1 
   vdt v_sum3a =  vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);//BODYCOM1 
   vdt v_sum4a =  vec_madd(v_xa, v_tauR_D2, v_tauR_D1);//BODYCOM1 
// 2//BODYCOM1 
   v_mhu_A1  = VEC_LDS(0, &mhu_a[4]);//BODYCOM1 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[5]);//BODYCOM1 
   v_mhu_B1 = VEC_LDS(0, &mhu_a[9]);//BODYCOM1 
   v_mhu_B2 = VEC_LDS(0, &mhu_a[10]);//BODYCOM1 
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[7]);//BODYCOM1 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[8]);//BODYCOM1 
  v_tauR_D1   = VEC_LDS(0, &tauR_a[20]);//BODYCOM1 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[21]);//BODYCOM1 
//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM1 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM1 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM1 
// 3//BODYCOM1 
//BODYCOM1 
   v_mhu_A1  = VEC_LDS(0, &mhu_a[2]);//BODYCOM1 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[3]);//BODYCOM1 
   v_mhu_B2 = VEC_LDS(0, &mhu_a[8]);//BODYCOM1 
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[5]);//BODYCOM1 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[6]);//BODYCOM1 
  v_tauR_D1   = VEC_LDS(0, &tauR_a[18]);//BODYCOM1 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[19]);//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM1 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM1 
// 4//BODYCOM1 
   v_mhu_A1  = VEC_LDS(0, &mhu_a[0]);//BODYCOM1 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[1]);//BODYCOM1 
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[3]);//BODYCOM1 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[4]);//BODYCOM1 
  v_tauR_D1   = VEC_LDS(0, &tauR_a[16]);//BODYCOM1 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[17]);//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM1 
// 5//BODYCOM1 
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[1]);//BODYCOM1 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[2]);//BODYCOM1 
  v_tauR_D1   = VEC_LDS(0, &tauR_a[14]);//BODYCOM1 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[15]);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM1 
// 6//BODYCOM1 
   v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[0]);//BODYCOM1 
   v_tauR_D1   = VEC_LDS(0, &tauR_a[12]);//BODYCOM1 
   v_tauR_D2   = VEC_LDS(0, &tauR_a[13]);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2); // Was missing//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM1 
// 7//BODYCOM1 
   v_tauR_D2   = VEC_LDS(0, &tauR_a[11]);//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM1 
//BODYCOM1 
// A//BODYCOM1 
   vdt temp1a   = vec_swdiv_nochk(v_sum3a,v_sum4a); //BODYCOM1 
   vdt temp2a   = vec_swdiv_nochk(v_sum1a,v_sum2a); //BODYCOM1 
   vdt temp3a   = vec_ld(0, &g[ii]);  //BODYCOM1 
   temp2a   = vec_sub(temp2a, temp3a);  //BODYCOM1 
   temp3a   = vec_madd(temp2a, temp1a, temp3a);//BODYCOM1 

   
   //vec_st(temp3a, 0, &g[ii]);
   switch (leftover) {
   case 1:
     g[ii+0] = temp3a get [0];
     break;

   case 2:
     g[ii+0] = temp3a get [0];
     g[ii+1] = temp3a get [1];
     break;

   case 3:
     g[ii+0] = temp3a get [0];
     g[ii+1] = temp3a get [1];
     g[ii+2] = temp3a get [2];
     break;

     default:
    //println("BAD LEFTOVER: %d\n", leftover); fflush(stdout);
    assert(0);
   }
 }

}



void update_f2Gate_v1(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a)
{
  //update_f2Gate(dt, nCells, VM, g, mhu_a, tauR_a); return;
  typedef vector4double vdt;


int gateIndex= 9;
int mhu_l= 8 ;
int mhu_m= 5;
int tauR_l=11;
int tauR_m=12;

/*
 int mhu_k  = 12;  
 int tauR_k = 22;  
*/
   
 double  tauRdt_a[12];
         tauRdt_a[0]  = tauR_a[0]*dt;
         tauRdt_a[1]  = tauR_a[1]*dt;
         tauRdt_a[2]  = tauR_a[2]*dt;
         tauRdt_a[3]  = tauR_a[3]*dt;
         tauRdt_a[4]  = tauR_a[4]*dt;
         tauRdt_a[5]  = tauR_a[5]*dt;
         tauRdt_a[6]  = tauR_a[6]*dt; 
         tauRdt_a[7]  = tauR_a[7]*dt;
         tauRdt_a[8]  = tauR_a[8]*dt;
         tauRdt_a[9]  = tauR_a[9]*dt;
         tauRdt_a[10] = tauR_a[10]*dt;
         tauRdt_a[11] = tauR_a[11]*dt;

 int ii=0;

 // 0 aligned
 // 1: load produces (JUNK, &VM[0], &VM[1], &VM[2])
 // 2: load produces (JUNK, JUNK, &VM[0], &VM[1])
 // 3: load produces (JUNK, JUNK, JUNK, &VM[0])
 int misalignment = ((long int)&VM[0] % 32) / 8;
 assert(misalignment == 0);

 if (misalignment) { 
      vdt v_xa   = vec_ld(0, &VM[ii]);  //BODYCOM0 
//BODYCOM0 
// 1//BODYCOM0 
 /*//BODYCOM0 
   sum1[ii] = mhu_a[4];//BODYCOM0 
   sum1[ii] = mhu_a[3]     + VM[ii]*sum1[ii];//BODYCOM0 
   sum2[ii] = mhu_a[12];//BODYCOM0 
   sum2[ii] = mhu_a[11]    + VM[ii]*sum2[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[11];//BODYCOM0 
   sum3[ii] = tauRdt_a[10] + VM[ii]*sum3[ii];//BODYCOM0 
   sum4[ii] = tauR_a[22]; //BODYCOM0 
   sum4[ii] = tauR_a[21]   + VM[ii]*sum4[ii]; //BODYCOM0 
*///BODYCOM0 
//BODYCOM0 
  vdt v_mhu_A1    = VEC_LDS(0, &mhu_a[3]);//BODYCOM0 
  vdt v_mhu_A2    = VEC_LDS(0, &mhu_a[4]);//BODYCOM0 
  vdt v_mhu_B1    = VEC_LDS(0, &mhu_a[11]);//BODYCOM0 
  vdt v_mhu_B2    = VEC_LDS(0, &mhu_a[12]);//BODYCOM0 
  vdt v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[10]);//BODYCOM0 
  vdt v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[11]);//BODYCOM0 
  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[21]);//BODYCOM0 
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[22]);//BODYCOM0 
//BODYCOM0 
  vdt v_sum1a = vec_madd(v_xa, v_mhu_A2,    v_mhu_A1);//BODYCOM0 
  vdt v_sum2a = vec_madd(v_xa, v_mhu_B2,    v_mhu_B1);//BODYCOM0 
  vdt v_sum3a = vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);//BODYCOM0 
  vdt v_sum4a = vec_madd(v_xa, v_tauR_D2,   v_tauR_D1);//BODYCOM0 
//BODYCOM0 
// 2  //BODYCOM0 
/*//BODYCOM0 
   sum1[ii] = mhu_a[2]    + VM[ii]*sum1[ii];//BODYCOM0 
   sum1[ii] = mhu_a[1]    + VM[ii]*sum1[ii];//BODYCOM0 
   sum2[ii] = mhu_a[10]   + VM[ii]*sum2[ii];//BODYCOM0 
   sum2[ii] = mhu_a[9]    + VM[ii]*sum2[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[9] + VM[ii]*sum3[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[8] + VM[ii]*sum3[ii];//BODYCOM0 
   sum4[ii] = tauR_a[20]  + VM[ii]*sum4[ii]; //BODYCOM0 
   sum4[ii] = tauR_a[19]  + VM[ii]*sum4[ii]; //BODYCOM0 
*///BODYCOM0 
   v_mhu_A1    = VEC_LDS(0, &mhu_a[1]);//BODYCOM0 
   v_mhu_A2    = VEC_LDS(0, &mhu_a[2]);//BODYCOM0 
   v_mhu_B1    = VEC_LDS(0, &mhu_a[9]);//BODYCOM0 
   v_mhu_B2    = VEC_LDS(0, &mhu_a[10]);//BODYCOM0 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[8]);//BODYCOM0 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[9]);//BODYCOM0 
   v_tauR_D1   = VEC_LDS(0, &tauR_a[19]);//BODYCOM0 
   v_tauR_D2   = VEC_LDS(0, &tauR_a[20]);//BODYCOM0 
//BODYCOM0 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM0 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM0 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM0 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM0 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM0 
//BODYCOM0 
// 3//BODYCOM0 
/*//BODYCOM0 
   sum1[ii] = mhu_a[0]    + VM[ii]*sum1[ii];//BODYCOM0 
   sum2[ii] = mhu_a[8]    + VM[ii]*sum2[ii];//BODYCOM0 
   sum2[ii] = mhu_a[7]    + VM[ii]*sum2[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[7] + VM[ii]*sum3[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[6] + VM[ii]*sum3[ii];//BODYCOM0 
   sum4[ii] = tauR_a[18]  + VM[ii]*sum4[ii]; //BODYCOM0 
   sum4[ii] = tauR_a[17]  + VM[ii]*sum4[ii]; //BODYCOM0 
*///BODYCOM0 
//BODYCOM0 
//BODYCOM0 
   v_mhu_A1    = VEC_LDS(0, &mhu_a[0]);//BODYCOM0 
//BODYCOM0 
   v_mhu_B1    = VEC_LDS(0, &mhu_a[7]);//BODYCOM0 
   v_mhu_B2    = VEC_LDS(0, &mhu_a[8]);//BODYCOM0 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[6]);//BODYCOM0 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[7]);//BODYCOM0 
   v_tauR_D1   = VEC_LDS(0, &tauR_a[17]);//BODYCOM0 
   v_tauR_D2   = VEC_LDS(0, &tauR_a[18]);//BODYCOM0 
//BODYCOM0 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM0 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM0 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM0 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM0 
//BODYCOM0 
// 4//BODYCOM0 
/*//BODYCOM0 
   sum2[ii] = mhu_a[6]    + VM[ii]*sum2[ii];//BODYCOM0 
   sum2[ii] = mhu_a[5]    + VM[ii]*sum2[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[5] + VM[ii]*sum3[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[4] + VM[ii]*sum3[ii];//BODYCOM0 
   sum4[ii] = tauR_a[16]  + VM[ii]*sum4[ii]; //BODYCOM0 
   sum4[ii] = tauR_a[15]  + VM[ii]*sum4[ii]; //BODYCOM0 
*///BODYCOM0 
//BODYCOM0 
   v_mhu_B1    = VEC_LDS(0, &mhu_a[5]);//BODYCOM0 
   v_mhu_B2    = VEC_LDS(0, &mhu_a[6]);//BODYCOM0 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[4]);//BODYCOM0 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[5]);//BODYCOM0 
   v_tauR_D1   = VEC_LDS(0, &tauR_a[15]);//BODYCOM0 
   v_tauR_D2   = VEC_LDS(0, &tauR_a[16]);//BODYCOM0 
//BODYCOM0 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM0 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM0 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM0 
//BODYCOM0 
// 5//BODYCOM0 
/* //BODYCOM0 
   sum3[ii] = tauRdt_a[3] + VM[ii]*sum3[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[2] + VM[ii]*sum3[ii];//BODYCOM0 
   sum4[ii] = tauR_a[14]  + VM[ii]*sum4[ii]; //BODYCOM0 
   sum4[ii] = tauR_a[13]  + VM[ii]*sum4[ii]; //BODYCOM0 
*///BODYCOM0 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[2]);//BODYCOM0 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[3]);//BODYCOM0 
   v_tauR_D1   = VEC_LDS(0, &tauR_a[13]);//BODYCOM0 
   v_tauR_D2   = VEC_LDS(0, &tauR_a[14]);//BODYCOM0 
//BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM0 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM0 
//BODYCOM0 
// 6//BODYCOM0 
/*//BODYCOM0 
   sum3[ii] = tauRdt_a[1] + VM[ii]*sum3[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[0] + VM[ii]*sum3[ii];//BODYCOM0 
   sum4[ii] = tauR_a[12]  + VM[ii]*sum4[ii]; //BODYCOM0 
*///BODYCOM0 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[0]);//BODYCOM0 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[1]);//BODYCOM0 
   v_tauR_D1   = VEC_LDS(0, &tauR_a[12]);//BODYCOM0 
//BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM0 
//BODYCOM0 
/*//BODYCOM0 
   tauRdt[ii] = sum3[ii]/sum4[ii];//BODYCOM0 
   mhu[ii]    = sum1[ii]/sum2[ii];//BODYCOM0 
//BODYCOM0 
//BODYCOM0 
//BODYCOM0 
   g[ii] +=  mhu[ii]*tauRdt[ii] - g[ii]*tauRdt[ii];//BODYCOM0 
*///BODYCOM0 
//BODYCOM0 
// Still need to reuse the register I have above ... the code will look like gibberish because//BODYCOM0 
// the targets won't make sense.//BODYCOM0 
// A//BODYCOM0 
   v_xa      = vec_swdiv_nochk(v_sum3a, v_sum4a); //BODYCOM0 
   v_mhu_A2  = vec_swdiv_nochk(v_sum1a, v_sum2a); //BODYCOM0 
   v_mhu_B2  = vec_ld(0, &g[ii]);  //BODYCOM0 
   v_tauR_D2 = vec_sub(v_mhu_A2, v_mhu_B2);  //BODYCOM0 
   v_mhu_B2  = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);//BODYCOM0 

   
   //vec_st(v_mhu_B2, 0, &g[ii]);
   switch (misalignment) {
   case 1:
     g[ii+0] = v_mhu_B2 get [1];
     g[ii+1] = v_mhu_B2 get [2];
     g[ii+2] = v_mhu_B2 get [3];
     break;

   case 2:
     g[ii+0] = v_mhu_B2 get [2];
     g[ii+1] = v_mhu_B2 get [3];
     break;

   case 3:
     g[ii+0] = v_mhu_B2 get [3];
     break;

     default:
    //println("BAD MISALIGNMENT: %d\n", misalignment); fflush(stdout);
     assert(0);
   }

   ii = 4 - misalignment;
 }

 for (ii;ii<TTROUND(nCells,8);ii+=8)
 { 
   vdt v_xa   = vec_ld(0, &VM[ii]);  
   vdt v_xb   = vec_ld(0, &VM[ii+4]);  

// 1
 /*
   sum1[ii] = mhu_a[4];
   sum1[ii] = mhu_a[3]     + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[12];
   sum2[ii] = mhu_a[11]    + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[11];
   sum3[ii] = tauRdt_a[10] + VM[ii]*sum3[ii];
   sum4[ii] = tauR_a[22]; 
   sum4[ii] = tauR_a[21]   + VM[ii]*sum4[ii]; 
*/

  vdt v_mhu_A1    = VEC_LDS(0, &mhu_a[3]);
  vdt v_mhu_A2    = VEC_LDS(0, &mhu_a[4]);
  vdt v_mhu_B1    = VEC_LDS(0, &mhu_a[11]);
  vdt v_mhu_B2    = VEC_LDS(0, &mhu_a[12]);
  vdt v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[10]);
  vdt v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[11]);
  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[21]);
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[22]);

  vdt v_sum1a = vec_madd(v_xa, v_mhu_A2,    v_mhu_A1);
  vdt v_sum1b = vec_madd(v_xb, v_mhu_A2,    v_mhu_A1);
  vdt v_sum2a = vec_madd(v_xa, v_mhu_B2,    v_mhu_B1);
  vdt v_sum2b = vec_madd(v_xb, v_mhu_B2,    v_mhu_B1);
  vdt v_sum3a = vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);
  vdt v_sum3b = vec_madd(v_xb, v_tauRdt_C2, v_tauRdt_C1);
  vdt v_sum4a = vec_madd(v_xa, v_tauR_D2,   v_tauR_D1);
  vdt v_sum4b = vec_madd(v_xb, v_tauR_D2,   v_tauR_D1);

// 2  
/*
   sum1[ii] = mhu_a[2]    + VM[ii]*sum1[ii];
   sum1[ii] = mhu_a[1]    + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[10]   + VM[ii]*sum2[ii];
   sum2[ii] = mhu_a[9]    + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[9] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[8] + VM[ii]*sum3[ii];
   sum4[ii] = tauR_a[20]  + VM[ii]*sum4[ii]; 
   sum4[ii] = tauR_a[19]  + VM[ii]*sum4[ii]; 
*/
   v_mhu_A1    = VEC_LDS(0, &mhu_a[1]);
   v_mhu_A2    = VEC_LDS(0, &mhu_a[2]);
   v_mhu_B1    = VEC_LDS(0, &mhu_a[9]);
   v_mhu_B2    = VEC_LDS(0, &mhu_a[10]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[8]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[9]);
   v_tauR_D1   = VEC_LDS(0, &tauR_a[19]);
   v_tauR_D2   = VEC_LDS(0, &tauR_a[20]);

   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);

   v_sum1b = vec_madd(v_xb, v_sum1b, v_mhu_A2);
   v_sum1b = vec_madd(v_xb, v_sum1b, v_mhu_A1);
   v_sum2b = vec_madd(v_xb, v_sum2b, v_mhu_B2);
   v_sum2b = vec_madd(v_xb, v_sum2b, v_mhu_B1); 
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   v_sum4b = vec_madd(v_xb, v_sum4b, v_tauR_D2);
   v_sum4b = vec_madd(v_xb, v_sum4b, v_tauR_D1);

// 3
/*
   sum1[ii] = mhu_a[0]    + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[8]    + VM[ii]*sum2[ii];
   sum2[ii] = mhu_a[7]    + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[7] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[6] + VM[ii]*sum3[ii];
   sum4[ii] = tauR_a[18]  + VM[ii]*sum4[ii]; 
   sum4[ii] = tauR_a[17]  + VM[ii]*sum4[ii]; 
*/


   v_mhu_A1    = VEC_LDS(0, &mhu_a[0]);

   v_mhu_B1    = VEC_LDS(0, &mhu_a[7]);
   v_mhu_B2    = VEC_LDS(0, &mhu_a[8]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[6]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[7]);
   v_tauR_D1   = VEC_LDS(0, &tauR_a[17]);
   v_tauR_D2   = VEC_LDS(0, &tauR_a[18]);

   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);

   v_sum1b = vec_madd(v_xb, v_sum1b, v_mhu_A1);
   v_sum2b = vec_madd(v_xb, v_sum2b, v_mhu_B2);
   v_sum2b = vec_madd(v_xb, v_sum2b, v_mhu_B1); 
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   v_sum4b = vec_madd(v_xb, v_sum4b, v_tauR_D2);
   v_sum4b = vec_madd(v_xb, v_sum4b, v_tauR_D1);

// 4
/*
   sum2[ii] = mhu_a[6]    + VM[ii]*sum2[ii];
   sum2[ii] = mhu_a[5]    + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[5] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[4] + VM[ii]*sum3[ii];
   sum4[ii] = tauR_a[16]  + VM[ii]*sum4[ii]; 
   sum4[ii] = tauR_a[15]  + VM[ii]*sum4[ii]; 
*/

   v_mhu_B1    = VEC_LDS(0, &mhu_a[5]);
   v_mhu_B2    = VEC_LDS(0, &mhu_a[6]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[4]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[5]);
   v_tauR_D1   = VEC_LDS(0, &tauR_a[15]);
   v_tauR_D2   = VEC_LDS(0, &tauR_a[16]);

   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);

   v_sum2b = vec_madd(v_xb, v_sum2b, v_mhu_B2);
   v_sum2b = vec_madd(v_xb, v_sum2b, v_mhu_B1); 
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   v_sum4b = vec_madd(v_xb, v_sum4b, v_tauR_D2);
   v_sum4b = vec_madd(v_xb, v_sum4b, v_tauR_D1);


// 5
/* 
   sum3[ii] = tauRdt_a[3] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[2] + VM[ii]*sum3[ii];
   sum4[ii] = tauR_a[14]  + VM[ii]*sum4[ii]; 
   sum4[ii] = tauR_a[13]  + VM[ii]*sum4[ii]; 
*/
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[2]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[3]);
   v_tauR_D1   = VEC_LDS(0, &tauR_a[13]);
   v_tauR_D2   = VEC_LDS(0, &tauR_a[14]);

   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);

   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   v_sum4b = vec_madd(v_xb, v_sum4b, v_tauR_D2);
   v_sum4b = vec_madd(v_xb, v_sum4b, v_tauR_D1);

// 6
/*
   sum3[ii] = tauRdt_a[1] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[0] + VM[ii]*sum3[ii];
   sum4[ii] = tauR_a[12]  + VM[ii]*sum4[ii]; 
*/
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[0]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[1]);
   v_tauR_D1   = VEC_LDS(0, &tauR_a[12]);

   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);

   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   v_sum4b = vec_madd(v_xb, v_sum4b, v_tauR_D1);

/*
   tauRdt[ii] = sum3[ii]/sum4[ii];
   mhu[ii]    = sum1[ii]/sum2[ii];



   g[ii] +=  mhu[ii]*tauRdt[ii] - g[ii]*tauRdt[ii];
*/

// Still need to reuse the register I have above ... the code will look like gibberish because
// the targets won't make sense.
// A
   v_xa      = vec_swdiv_nochk(v_sum3a, v_sum4a); 
   v_mhu_A2  = vec_swdiv_nochk(v_sum1a, v_sum2a); 
   v_mhu_B2  = vec_ld(0, &g[ii]);  
   v_tauR_D2 = vec_sub(v_mhu_A2, v_mhu_B2);  
   v_mhu_B2  = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);
   vec_st(v_mhu_B2, 0, &g[ii]);

// B
   v_xb      = vec_swdiv_nochk(v_sum3b, v_sum4b); 
   v_mhu_A1  = vec_swdiv_nochk(v_sum1b, v_sum2b); 
   v_mhu_B1  = vec_ld(0, &g[ii+4]);  
   v_tauR_D1 = vec_sub(v_mhu_A1, v_mhu_B1);  
   v_mhu_B1  = vec_madd(v_tauR_D1, v_xb, v_mhu_B1);
   vec_st(v_mhu_B1, 0, &g[ii+4]);
 } 

 for (ii;ii<TTROUND(nCells,4);ii+=4)
 { 
 //BODYSTART
   vdt v_xa   = vec_ld(0, &VM[ii]);  

// 1
 /*
   sum1[ii] = mhu_a[4];
   sum1[ii] = mhu_a[3]     + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[12];
   sum2[ii] = mhu_a[11]    + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[11];
   sum3[ii] = tauRdt_a[10] + VM[ii]*sum3[ii];
   sum4[ii] = tauR_a[22]; 
   sum4[ii] = tauR_a[21]   + VM[ii]*sum4[ii]; 
*/

  vdt v_mhu_A1    = VEC_LDS(0, &mhu_a[3]);
  vdt v_mhu_A2    = VEC_LDS(0, &mhu_a[4]);
  vdt v_mhu_B1    = VEC_LDS(0, &mhu_a[11]);
  vdt v_mhu_B2    = VEC_LDS(0, &mhu_a[12]);
  vdt v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[10]);
  vdt v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[11]);
  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[21]);
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[22]);

  vdt v_sum1a = vec_madd(v_xa, v_mhu_A2,    v_mhu_A1);
  vdt v_sum2a = vec_madd(v_xa, v_mhu_B2,    v_mhu_B1);
  vdt v_sum3a = vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);
  vdt v_sum4a = vec_madd(v_xa, v_tauR_D2,   v_tauR_D1);

// 2  
/*
   sum1[ii] = mhu_a[2]    + VM[ii]*sum1[ii];
   sum1[ii] = mhu_a[1]    + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[10]   + VM[ii]*sum2[ii];
   sum2[ii] = mhu_a[9]    + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[9] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[8] + VM[ii]*sum3[ii];
   sum4[ii] = tauR_a[20]  + VM[ii]*sum4[ii]; 
   sum4[ii] = tauR_a[19]  + VM[ii]*sum4[ii]; 
*/
   v_mhu_A1    = VEC_LDS(0, &mhu_a[1]);
   v_mhu_A2    = VEC_LDS(0, &mhu_a[2]);
   v_mhu_B1    = VEC_LDS(0, &mhu_a[9]);
   v_mhu_B2    = VEC_LDS(0, &mhu_a[10]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[8]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[9]);
   v_tauR_D1   = VEC_LDS(0, &tauR_a[19]);
   v_tauR_D2   = VEC_LDS(0, &tauR_a[20]);

   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);

// 3
/*
   sum1[ii] = mhu_a[0]    + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[8]    + VM[ii]*sum2[ii];
   sum2[ii] = mhu_a[7]    + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[7] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[6] + VM[ii]*sum3[ii];
   sum4[ii] = tauR_a[18]  + VM[ii]*sum4[ii]; 
   sum4[ii] = tauR_a[17]  + VM[ii]*sum4[ii]; 
*/


   v_mhu_A1    = VEC_LDS(0, &mhu_a[0]);

   v_mhu_B1    = VEC_LDS(0, &mhu_a[7]);
   v_mhu_B2    = VEC_LDS(0, &mhu_a[8]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[6]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[7]);
   v_tauR_D1   = VEC_LDS(0, &tauR_a[17]);
   v_tauR_D2   = VEC_LDS(0, &tauR_a[18]);

   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);

// 4
/*
   sum2[ii] = mhu_a[6]    + VM[ii]*sum2[ii];
   sum2[ii] = mhu_a[5]    + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[5] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[4] + VM[ii]*sum3[ii];
   sum4[ii] = tauR_a[16]  + VM[ii]*sum4[ii]; 
   sum4[ii] = tauR_a[15]  + VM[ii]*sum4[ii]; 
*/

   v_mhu_B1    = VEC_LDS(0, &mhu_a[5]);
   v_mhu_B2    = VEC_LDS(0, &mhu_a[6]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[4]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[5]);
   v_tauR_D1   = VEC_LDS(0, &tauR_a[15]);
   v_tauR_D2   = VEC_LDS(0, &tauR_a[16]);

   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);

// 5
/* 
   sum3[ii] = tauRdt_a[3] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[2] + VM[ii]*sum3[ii];
   sum4[ii] = tauR_a[14]  + VM[ii]*sum4[ii]; 
   sum4[ii] = tauR_a[13]  + VM[ii]*sum4[ii]; 
*/
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[2]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[3]);
   v_tauR_D1   = VEC_LDS(0, &tauR_a[13]);
   v_tauR_D2   = VEC_LDS(0, &tauR_a[14]);

   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);

// 6
/*
   sum3[ii] = tauRdt_a[1] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[0] + VM[ii]*sum3[ii];
   sum4[ii] = tauR_a[12]  + VM[ii]*sum4[ii]; 
*/
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[0]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[1]);
   v_tauR_D1   = VEC_LDS(0, &tauR_a[12]);

   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);

/*
   tauRdt[ii] = sum3[ii]/sum4[ii];
   mhu[ii]    = sum1[ii]/sum2[ii];



   g[ii] +=  mhu[ii]*tauRdt[ii] - g[ii]*tauRdt[ii];
*/

// Still need to reuse the register I have above ... the code will look like gibberish because
// the targets won't make sense.
// A
   v_xa      = vec_swdiv_nochk(v_sum3a, v_sum4a); 
   v_mhu_A2  = vec_swdiv_nochk(v_sum1a, v_sum2a); 
   v_mhu_B2  = vec_ld(0, &g[ii]);  
   v_tauR_D2 = vec_sub(v_mhu_A2, v_mhu_B2);  
   v_mhu_B2  = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);
   vec_st(v_mhu_B2, 0, &g[ii]); //BODYEND
 } 

 // the elements in the end of the vector should be ignored
 int leftover = nCells - ii;
 if (leftover) { 
      vdt v_xa   = vec_ld(0, &VM[ii]);  //BODYCOM1 
//BODYCOM1 
// 1//BODYCOM1 
 /*//BODYCOM1 
   sum1[ii] = mhu_a[4];//BODYCOM1 
   sum1[ii] = mhu_a[3]     + VM[ii]*sum1[ii];//BODYCOM1 
   sum2[ii] = mhu_a[12];//BODYCOM1 
   sum2[ii] = mhu_a[11]    + VM[ii]*sum2[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[11];//BODYCOM1 
   sum3[ii] = tauRdt_a[10] + VM[ii]*sum3[ii];//BODYCOM1 
   sum4[ii] = tauR_a[22]; //BODYCOM1 
   sum4[ii] = tauR_a[21]   + VM[ii]*sum4[ii]; //BODYCOM1 
*///BODYCOM1 
//BODYCOM1 
  vdt v_mhu_A1    = VEC_LDS(0, &mhu_a[3]);//BODYCOM1 
  vdt v_mhu_A2    = VEC_LDS(0, &mhu_a[4]);//BODYCOM1 
  vdt v_mhu_B1    = VEC_LDS(0, &mhu_a[11]);//BODYCOM1 
  vdt v_mhu_B2    = VEC_LDS(0, &mhu_a[12]);//BODYCOM1 
  vdt v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[10]);//BODYCOM1 
  vdt v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[11]);//BODYCOM1 
  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[21]);//BODYCOM1 
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[22]);//BODYCOM1 
//BODYCOM1 
  vdt v_sum1a = vec_madd(v_xa, v_mhu_A2,    v_mhu_A1);//BODYCOM1 
  vdt v_sum2a = vec_madd(v_xa, v_mhu_B2,    v_mhu_B1);//BODYCOM1 
  vdt v_sum3a = vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);//BODYCOM1 
  vdt v_sum4a = vec_madd(v_xa, v_tauR_D2,   v_tauR_D1);//BODYCOM1 
//BODYCOM1 
// 2  //BODYCOM1 
/*//BODYCOM1 
   sum1[ii] = mhu_a[2]    + VM[ii]*sum1[ii];//BODYCOM1 
   sum1[ii] = mhu_a[1]    + VM[ii]*sum1[ii];//BODYCOM1 
   sum2[ii] = mhu_a[10]   + VM[ii]*sum2[ii];//BODYCOM1 
   sum2[ii] = mhu_a[9]    + VM[ii]*sum2[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[9] + VM[ii]*sum3[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[8] + VM[ii]*sum3[ii];//BODYCOM1 
   sum4[ii] = tauR_a[20]  + VM[ii]*sum4[ii]; //BODYCOM1 
   sum4[ii] = tauR_a[19]  + VM[ii]*sum4[ii]; //BODYCOM1 
*///BODYCOM1 
   v_mhu_A1    = VEC_LDS(0, &mhu_a[1]);//BODYCOM1 
   v_mhu_A2    = VEC_LDS(0, &mhu_a[2]);//BODYCOM1 
   v_mhu_B1    = VEC_LDS(0, &mhu_a[9]);//BODYCOM1 
   v_mhu_B2    = VEC_LDS(0, &mhu_a[10]);//BODYCOM1 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[8]);//BODYCOM1 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[9]);//BODYCOM1 
   v_tauR_D1   = VEC_LDS(0, &tauR_a[19]);//BODYCOM1 
   v_tauR_D2   = VEC_LDS(0, &tauR_a[20]);//BODYCOM1 
//BODYCOM1 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM1 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM1 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM1 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM1 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM1 
//BODYCOM1 
// 3//BODYCOM1 
/*//BODYCOM1 
   sum1[ii] = mhu_a[0]    + VM[ii]*sum1[ii];//BODYCOM1 
   sum2[ii] = mhu_a[8]    + VM[ii]*sum2[ii];//BODYCOM1 
   sum2[ii] = mhu_a[7]    + VM[ii]*sum2[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[7] + VM[ii]*sum3[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[6] + VM[ii]*sum3[ii];//BODYCOM1 
   sum4[ii] = tauR_a[18]  + VM[ii]*sum4[ii]; //BODYCOM1 
   sum4[ii] = tauR_a[17]  + VM[ii]*sum4[ii]; //BODYCOM1 
*///BODYCOM1 
//BODYCOM1 
//BODYCOM1 
   v_mhu_A1    = VEC_LDS(0, &mhu_a[0]);//BODYCOM1 
//BODYCOM1 
   v_mhu_B1    = VEC_LDS(0, &mhu_a[7]);//BODYCOM1 
   v_mhu_B2    = VEC_LDS(0, &mhu_a[8]);//BODYCOM1 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[6]);//BODYCOM1 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[7]);//BODYCOM1 
   v_tauR_D1   = VEC_LDS(0, &tauR_a[17]);//BODYCOM1 
   v_tauR_D2   = VEC_LDS(0, &tauR_a[18]);//BODYCOM1 
//BODYCOM1 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM1 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM1 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM1 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM1 
//BODYCOM1 
// 4//BODYCOM1 
/*//BODYCOM1 
   sum2[ii] = mhu_a[6]    + VM[ii]*sum2[ii];//BODYCOM1 
   sum2[ii] = mhu_a[5]    + VM[ii]*sum2[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[5] + VM[ii]*sum3[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[4] + VM[ii]*sum3[ii];//BODYCOM1 
   sum4[ii] = tauR_a[16]  + VM[ii]*sum4[ii]; //BODYCOM1 
   sum4[ii] = tauR_a[15]  + VM[ii]*sum4[ii]; //BODYCOM1 
*///BODYCOM1 
//BODYCOM1 
   v_mhu_B1    = VEC_LDS(0, &mhu_a[5]);//BODYCOM1 
   v_mhu_B2    = VEC_LDS(0, &mhu_a[6]);//BODYCOM1 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[4]);//BODYCOM1 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[5]);//BODYCOM1 
   v_tauR_D1   = VEC_LDS(0, &tauR_a[15]);//BODYCOM1 
   v_tauR_D2   = VEC_LDS(0, &tauR_a[16]);//BODYCOM1 
//BODYCOM1 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM1 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM1 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM1 
//BODYCOM1 
// 5//BODYCOM1 
/* //BODYCOM1 
   sum3[ii] = tauRdt_a[3] + VM[ii]*sum3[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[2] + VM[ii]*sum3[ii];//BODYCOM1 
   sum4[ii] = tauR_a[14]  + VM[ii]*sum4[ii]; //BODYCOM1 
   sum4[ii] = tauR_a[13]  + VM[ii]*sum4[ii]; //BODYCOM1 
*///BODYCOM1 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[2]);//BODYCOM1 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[3]);//BODYCOM1 
   v_tauR_D1   = VEC_LDS(0, &tauR_a[13]);//BODYCOM1 
   v_tauR_D2   = VEC_LDS(0, &tauR_a[14]);//BODYCOM1 
//BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM1 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM1 
//BODYCOM1 
// 6//BODYCOM1 
/*//BODYCOM1 
   sum3[ii] = tauRdt_a[1] + VM[ii]*sum3[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[0] + VM[ii]*sum3[ii];//BODYCOM1 
   sum4[ii] = tauR_a[12]  + VM[ii]*sum4[ii]; //BODYCOM1 
*///BODYCOM1 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[0]);//BODYCOM1 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[1]);//BODYCOM1 
   v_tauR_D1   = VEC_LDS(0, &tauR_a[12]);//BODYCOM1 
//BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM1 
//BODYCOM1 
/*//BODYCOM1 
   tauRdt[ii] = sum3[ii]/sum4[ii];//BODYCOM1 
   mhu[ii]    = sum1[ii]/sum2[ii];//BODYCOM1 
//BODYCOM1 
//BODYCOM1 
//BODYCOM1 
   g[ii] +=  mhu[ii]*tauRdt[ii] - g[ii]*tauRdt[ii];//BODYCOM1 
*///BODYCOM1 
//BODYCOM1 
// Still need to reuse the register I have above ... the code will look like gibberish because//BODYCOM1 
// the targets won't make sense.//BODYCOM1 
// A//BODYCOM1 
   v_xa      = vec_swdiv_nochk(v_sum3a, v_sum4a); //BODYCOM1 
   v_mhu_A2  = vec_swdiv_nochk(v_sum1a, v_sum2a); //BODYCOM1 
   v_mhu_B2  = vec_ld(0, &g[ii]);  //BODYCOM1 
   v_tauR_D2 = vec_sub(v_mhu_A2, v_mhu_B2);  //BODYCOM1 
   v_mhu_B2  = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);//BODYCOM1 

   
   //vec_st(v_mhu_B2, 0, &g[ii]);
   //printf("%d: leftover=%d\n",gateIndex,leftover); 
   switch (leftover) {
   case 1:
     g[ii+0] = v_mhu_B2 get [0];
     break;

   case 2:
     g[ii+0] = v_mhu_B2 get [0];
     g[ii+1] = v_mhu_B2 get [1];
     break;

   case 3:
     g[ii+0] = v_mhu_B2 get [0];
     g[ii+1] = v_mhu_B2 get [1];
     g[ii+2] = v_mhu_B2 get [2];
     break;

     default:
    //println("BAD LEFTOVER: %d\n", leftover); fflush(stdout);
    assert(0);
   }
 }

  
}

void update_jLGate_v1(double dt, int nCells,  double *VM, double *g, double *mhu_a, double *tauR_a)
{
  //update_jLGate(dt, nCells,  VM, g, mhu_a, tauR_a); return;
  typedef vector4double vdt;

int gateIndex=10;
int mhu_l=6 ;
int mhu_m= 12;
int tauR_l=1;
int tauR_m=1;
   
  vdt v_mhu_a0  = VEC_LDS(0, &mhu_a[0]);
  vdt v_mhu_a1  = VEC_LDS(0, &mhu_a[1]);
  vdt v_mhu_a2  = VEC_LDS(0, &mhu_a[2]);
  vdt v_mhu_a3  = VEC_LDS(0, &mhu_a[3]);
  vdt v_mhu_a4  = VEC_LDS(0, &mhu_a[4]);
  vdt v_mhu_a5  = VEC_LDS(0, &mhu_a[5]);
  vdt v_mhu_a6  = VEC_LDS(0, &mhu_a[6]);
  vdt v_mhu_a7  = VEC_LDS(0, &mhu_a[7]);
  vdt v_mhu_a8  = VEC_LDS(0, &mhu_a[8]);
  vdt v_mhu_a9  = VEC_LDS(0, &mhu_a[9]);
  vdt v_mhu_a10  = VEC_LDS(0, &mhu_a[10]);
  vdt v_mhu_a11  = VEC_LDS(0, &mhu_a[11]);
  vdt v_mhu_a12  = VEC_LDS(0, &mhu_a[12]);
  vdt v_mhu_a13  = VEC_LDS(0, &mhu_a[13]);
  vdt v_mhu_a14  = VEC_LDS(0, &mhu_a[14]);
  vdt v_mhu_a15  = VEC_LDS(0, &mhu_a[15]);
  vdt v_mhu_a16  = VEC_LDS(0, &mhu_a[16]);
  vdt v_mhu_a17  = VEC_LDS(0, &mhu_a[17]);


  vdt temp1   = vec_splats(tauR_a[0]*dt);

  // vdt v_TWO = vec_splats(2.0);

 int mhu_k  = 17;  
 int tauR_k = 1;  
 int ii=0;

 // 0 aligned
 // 1: load produces (JUNK, &VM[0], &VM[1], &VM[2])
 // 2: load produces (JUNK, JUNK, &VM[0], &VM[1])
 // 3: load produces (JUNK, JUNK, JUNK, &VM[0])
 int misalignment = ((long int)&VM[0] % 32) / 8;
 assert(misalignment == 0);

 if (misalignment) { 
      vdt v_v = vec_ld(0, &VM[ii]); //BODYCOM0 
//BODYCOM0 
   vdt v_sum1   =  vec_madd(v_v, v_mhu_a11, v_mhu_a10);//BODYCOM0 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a9);//BODYCOM0 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a8);//BODYCOM0 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a7);//BODYCOM0 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a6);//BODYCOM0 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a5);//BODYCOM0 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a4); //BODYCOM0 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a3);//BODYCOM0 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a2);//BODYCOM0 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a1);//BODYCOM0 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a0);//BODYCOM0 
//BODYCOM0 
   vdt v_sum2   =  vec_madd(v_v, v_mhu_a17, v_mhu_a16);//BODYCOM0 
   v_sum2   =  vec_madd(v_v, v_sum2, v_mhu_a15);//BODYCOM0 
   v_sum2   =  vec_madd(v_v, v_sum2, v_mhu_a14); //BODYCOM0 
   v_sum2   =  vec_madd(v_v, v_sum2, v_mhu_a13);//BODYCOM0 
   v_sum2   =  vec_madd(v_v, v_sum2, v_mhu_a12);//BODYCOM0 
// A//BODYCOM0 
   // vdt temp1a   = vec_swdiv_nochk(v_sum3a,v_sum4a);  ... temp1 is as we want it//BODYCOM0 
   vdt temp2a   = vec_swdiv_nochk(v_sum1,v_sum2); //BODYCOM0 
   vdt temp3a   = vec_ld(0, &g[ii]);  //BODYCOM0 
   temp2a   = vec_sub(temp2a, temp3a);  //BODYCOM0 
   temp3a   = vec_madd(temp2a, temp1, temp3a);//BODYCOM0 

   
   //vec_st(temp3a, 0, &g[ii]);
   switch (misalignment) {
   case 1:
     g[ii+0] = temp3a get [1];
     g[ii+1] = temp3a get [2];
     g[ii+2] = temp3a get [3];
     break;

   case 2:
     g[ii+0] = temp3a get [2];
     g[ii+1] = temp3a get [3];
     break;

   case 3:
     g[ii+0] = temp3a get [3];
     break;

     default:
    //println("BAD MISALIGNMENT: %d\n", misalignment); fflush(stdout);
     assert(0);
   }

   ii = 4 - misalignment;
 }

 for (ii;ii<TTROUND(nCells,8);ii+=8)
 { 
   vdt v_v = vec_ld(0, &VM[ii]); 
   vdt v_v_4 = vec_ld(0, &VM[ii+4]);    

   vdt v_sum1   =  vec_madd(v_v, v_mhu_a11, v_mhu_a10);
   vdt v_sum1_4 =  vec_madd(v_v_4, v_mhu_a11, v_mhu_a10);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a9);
   v_sum1_4 =  vec_madd(v_v_4, v_sum1_4, v_mhu_a9);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a8);
   v_sum1_4 =  vec_madd(v_v_4, v_sum1_4, v_mhu_a8);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a7);
   v_sum1_4 =  vec_madd(v_v_4, v_sum1_4, v_mhu_a7);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a6);
   v_sum1_4 =  vec_madd(v_v_4, v_sum1_4, v_mhu_a6);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a5);
   v_sum1_4 =  vec_madd(v_v_4, v_sum1_4, v_mhu_a5);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a4); 
   v_sum1_4 =  vec_madd(v_v_4, v_sum1_4, v_mhu_a4); 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a3);
   v_sum1_4 =  vec_madd(v_v_4, v_sum1_4, v_mhu_a3);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a2);
   v_sum1_4 =  vec_madd(v_v_4, v_sum1_4, v_mhu_a2);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a1);
   v_sum1_4 =  vec_madd(v_v_4, v_sum1_4, v_mhu_a1);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a0);
   v_sum1_4 =  vec_madd(v_v_4, v_sum1_4, v_mhu_a0);

   vdt v_sum2   =  vec_madd(v_v, v_mhu_a17, v_mhu_a16);
   vdt v_sum2_4 =  vec_madd(v_v_4, v_mhu_a17, v_mhu_a16);
   v_sum2   =  vec_madd(v_v, v_sum2, v_mhu_a15);
   v_sum2_4 =  vec_madd(v_v_4, v_sum2_4, v_mhu_a15);
   v_sum2   =  vec_madd(v_v, v_sum2, v_mhu_a14); 
   v_sum2_4 =  vec_madd(v_v_4, v_sum2_4, v_mhu_a14); 
   v_sum2   =  vec_madd(v_v, v_sum2, v_mhu_a13);
   v_sum2_4 =  vec_madd(v_v_4, v_sum2_4, v_mhu_a13);
   v_sum2   =  vec_madd(v_v, v_sum2, v_mhu_a12);
   v_sum2_4 =  vec_madd(v_v_4, v_sum2_4, v_mhu_a12);
// A
   // vdt temp1a   = vec_swdiv_nochk(v_sum3a,v_sum4a);  ... temp1 is as we want it
   vdt temp2a   = vec_swdiv_nochk(v_sum1,v_sum2); 
   vdt temp3a   = vec_ld(0, &g[ii]);  
   temp2a   = vec_sub(temp2a, temp3a);  
   temp3a   = vec_madd(temp2a, temp1, temp3a);
   vec_st(temp3a, 0, &g[ii]);
   /* This can also be implemented as the following, which has a higher flop rate and is more straightforward
   // temp1  = vec_swdiv_nochk(v_sum3a,v_sum4a);  ... temp1 is as we want it
   temp2a  = vec_swdiv_nochk(v_sum1a,v_sum2a); 
   temp3a  = vec_ld(0, &g[ii]);
   temp2a   = vec_madd(temp1, temp2a, temp3a);
   temp3a  = vec_nmadd(temp1, temp3a, temp2a );
   vec_st(temp3a, 0, &g[ii]); */
// B
   // vdt temp1b   = vec_swdiv_nochk(v_sum3b,v_sum4b);  ... temp1 is as we want it
   vdt temp2b   = vec_swdiv_nochk(v_sum1_4,v_sum2_4); 
   vdt temp3b   = vec_ld(0, &g[ii+4]);  
   temp2b   = vec_sub(temp2b, temp3b);  
   temp3b   = vec_madd(temp2b, temp1, temp3b);
   vec_st(temp3b, 0, &g[ii+4]);
 }
 for (ii;ii<TTROUND(nCells,4);ii+=4)
 { 
 //BODYSTART
   vdt v_v = vec_ld(0, &VM[ii]); 

   vdt v_sum1   =  vec_madd(v_v, v_mhu_a11, v_mhu_a10);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a9);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a8);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a7);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a6);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a5);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a4); 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a3);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a2);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a1);
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a0);

   vdt v_sum2   =  vec_madd(v_v, v_mhu_a17, v_mhu_a16);
   v_sum2   =  vec_madd(v_v, v_sum2, v_mhu_a15);
   v_sum2   =  vec_madd(v_v, v_sum2, v_mhu_a14); 
   v_sum2   =  vec_madd(v_v, v_sum2, v_mhu_a13);
   v_sum2   =  vec_madd(v_v, v_sum2, v_mhu_a12);
// A
   // vdt temp1a   = vec_swdiv_nochk(v_sum3a,v_sum4a);  ... temp1 is as we want it
   vdt temp2a   = vec_swdiv_nochk(v_sum1,v_sum2); 
   vdt temp3a   = vec_ld(0, &g[ii]);  
   temp2a   = vec_sub(temp2a, temp3a);  
   temp3a   = vec_madd(temp2a, temp1, temp3a);
   vec_st(temp3a, 0, &g[ii]); //BODYEND
   /* This can also be implemented as the following, which has a higher flop rate and is more straightforward
   // temp1  = vec_swdiv_nochk(v_sum3a,v_sum4a);  ... temp1 is as we want it
   temp2a  = vec_swdiv_nochk(v_sum1a,v_sum2a); 
   temp3a  = vec_ld(0, &g[ii]);
   temp2a   = vec_madd(temp1, temp2a, temp3a);
   temp3a  = vec_nmadd(temp1, temp3a, temp2a );
   vec_st(temp3a, 0, &g[ii]); */
 }

 // the elements in the end of the vector should be ignored
 int leftover = nCells - ii;
 if (leftover) { 
      vdt v_v = vec_ld(0, &VM[ii]); //BODYCOM1 
//BODYCOM1 
   vdt v_sum1   =  vec_madd(v_v, v_mhu_a11, v_mhu_a10);//BODYCOM1 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a9);//BODYCOM1 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a8);//BODYCOM1 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a7);//BODYCOM1 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a6);//BODYCOM1 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a5);//BODYCOM1 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a4); //BODYCOM1 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a3);//BODYCOM1 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a2);//BODYCOM1 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a1);//BODYCOM1 
   v_sum1   =  vec_madd(v_v, v_sum1, v_mhu_a0);//BODYCOM1 
//BODYCOM1 
   vdt v_sum2   =  vec_madd(v_v, v_mhu_a17, v_mhu_a16);//BODYCOM1 
   v_sum2   =  vec_madd(v_v, v_sum2, v_mhu_a15);//BODYCOM1 
   v_sum2   =  vec_madd(v_v, v_sum2, v_mhu_a14); //BODYCOM1 
   v_sum2   =  vec_madd(v_v, v_sum2, v_mhu_a13);//BODYCOM1 
   v_sum2   =  vec_madd(v_v, v_sum2, v_mhu_a12);//BODYCOM1 
// A//BODYCOM1 
   // vdt temp1a   = vec_swdiv_nochk(v_sum3a,v_sum4a);  ... temp1 is as we want it//BODYCOM1 
   vdt temp2a   = vec_swdiv_nochk(v_sum1,v_sum2); //BODYCOM1 
   vdt temp3a   = vec_ld(0, &g[ii]);  //BODYCOM1 
   temp2a   = vec_sub(temp2a, temp3a);  //BODYCOM1 
   temp3a   = vec_madd(temp2a, temp1, temp3a);//BODYCOM1 

   //printf("%d: leftover=%d\n",gateIndex,leftover); 
   
   //vec_st(temp3a, 0, &g[ii]);
   switch (leftover) {
   case 1:
     g[ii+0] = temp3a get [0];
     break;

   case 2:
     g[ii+0] = temp3a get [0];
     g[ii+1] = temp3a get [1];
     break;

   case 3:
     g[ii+0] = temp3a get [0];
     g[ii+1] = temp3a get [1];
     g[ii+2] = temp3a get [2];
     break;

     default:
    //println("BAD LEFTOVER: %d\n", leftover); fflush(stdout);
    assert(0);
   }
 }

}



void update_s0Gate_v1(double dt, int nCells,  double *VM, double *g, double *mhu_a, double *tauR_a)
{
  //update_s0Gate(dt, nCells,  VM, g, mhu_a, tauR_a); return;
  typedef vector4double vdt;

int gateIndex= 11;
int mhu_l= 7 ;
int mhu_m= 8;
int tauR_l=6 ;
int tauR_m=7 ;

/*
 int mhu_k  = 14;  
 int tauR_k = 12;  
*/
   
 double  tauRdt_a[7];
// for (int j=6;j>=0;j--)    tauRdt_a[j] = tauR_a[j]*dt;
         tauRdt_a[0]  = tauR_a[0]*dt,
         tauRdt_a[1]  = tauR_a[1]*dt,
         tauRdt_a[2]  = tauR_a[2]*dt,
         tauRdt_a[3]  = tauR_a[3]*dt,
         tauRdt_a[4]  = tauR_a[4]*dt,
         tauRdt_a[5]  = tauR_a[5]*dt,
         tauRdt_a[6]  = tauR_a[6]*dt; 

 int ii=0;

 // 0 aligned
 // 1: load produces (JUNK, &VM[0], &VM[1], &VM[2])
 // 2: load produces (JUNK, JUNK, &VM[0], &VM[1])
 // 3: load produces (JUNK, JUNK, JUNK, &VM[0])
 int misalignment = ((long int)&VM[0] % 32) / 8;
 assert(misalignment == 0);

 if (misalignment) { 
      vdt v_xa   = vec_ld(0, &VM[ii]);  //BODYCOM0 
//BODYCOM0 
//1//BODYCOM0 
/*//BODYCOM0 
//1//BODYCOM0 
   sum1[ii] = mhu_a[7];//BODYCOM0 
   sum1[ii] = mhu_a[6]    + VM[ii]*sum1[ii];//BODYCOM0 
   sum2[ii] = mhu_a[14];//BODYCOM0 
   sum2[ii] = mhu_a[13]   + VM[ii]*sum2[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[6];//BODYCOM0 
   sum3[ii] = tauRdt_a[5] + VM[ii]*sum3[ii];//BODYCOM0 
   sum4[ii] = tauR_a[12]; //BODYCOM0 
   sum4[ii] = tauR_a[11]  + VM[ii]*sum4[ii]; //BODYCOM0 
*///BODYCOM0 
//BODYCOM0 
  vdt v_mhu_A1    = VEC_LDS(0, &mhu_a[6]);//BODYCOM0 
  vdt v_mhu_A2    = VEC_LDS(0, &mhu_a[7]);//BODYCOM0 
  vdt v_mhu_B1    = VEC_LDS(0, &mhu_a[13]);//BODYCOM0 
  vdt v_mhu_B2    = VEC_LDS(0, &mhu_a[14]);//BODYCOM0 
  vdt v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[5]);//BODYCOM0 
  vdt v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[6]);//BODYCOM0 
  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[11]);//BODYCOM0 
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[12]);//BODYCOM0 
//BODYCOM0 
  vdt v_sum1a = vec_madd(v_xa, v_mhu_A2,    v_mhu_A1);//BODYCOM0 
  vdt v_sum2a = vec_madd(v_xa, v_mhu_B2,    v_mhu_B1);//BODYCOM0 
  vdt v_sum3a = vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);//BODYCOM0 
  vdt v_sum4a = vec_madd(v_xa, v_tauR_D2,   v_tauR_D1);//BODYCOM0 
//BODYCOM0 
//2//BODYCOM0 
/*//BODYCOM0 
//BODYCOM0 
   sum1[ii] = mhu_a[5]    + VM[ii]*sum1[ii];//BODYCOM0 
   sum1[ii] = mhu_a[4]    + VM[ii]*sum1[ii];//BODYCOM0 
   sum2[ii] = mhu_a[12]   + VM[ii]*sum2[ii];//BODYCOM0 
   sum2[ii] = mhu_a[11]   + VM[ii]*sum2[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[4] + VM[ii]*sum3[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[3] + VM[ii]*sum3[ii];//BODYCOM0 
   sum4[ii] = tauR_a[10]  + VM[ii]*sum4[ii]; //BODYCOM0 
   sum4[ii] = tauR_a[9]   + VM[ii]*sum4[ii]; //BODYCOM0 
*///BODYCOM0 
//BODYCOM0 
   v_mhu_A1    = VEC_LDS(0, &mhu_a[4]);//BODYCOM0 
   v_mhu_A2    = VEC_LDS(0, &mhu_a[5]);//BODYCOM0 
   v_mhu_B1    = VEC_LDS(0, &mhu_a[11]);//BODYCOM0 
   v_mhu_B2    = VEC_LDS(0, &mhu_a[12]);//BODYCOM0 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[3]);//BODYCOM0 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[4]);//BODYCOM0 
   v_tauR_D1   = VEC_LDS(0, &tauR_a[9]);//BODYCOM0 
   v_tauR_D2   = VEC_LDS(0, &tauR_a[10]);//BODYCOM0 
//BODYCOM0 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM0 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM0 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM0 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM0 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM0 
//BODYCOM0 
//3//BODYCOM0 
/*//BODYCOM0 
   sum1[ii] = mhu_a[3]    + VM[ii]*sum1[ii];//BODYCOM0 
   sum1[ii] = mhu_a[2]    + VM[ii]*sum1[ii];//BODYCOM0 
   sum2[ii] = mhu_a[10]   + VM[ii]*sum2[ii];//BODYCOM0 
   sum2[ii] = mhu_a[9]    + VM[ii]*sum2[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[2] + VM[ii]*sum3[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[1] + VM[ii]*sum3[ii];//BODYCOM0 
   sum4[ii] = tauR_a[8]   + VM[ii]*sum4[ii]; //BODYCOM0 
   sum4[ii] = tauR_a[7]   + VM[ii]*sum4[ii]; //BODYCOM0 
*///BODYCOM0 
//BODYCOM0 
   v_mhu_A1    = VEC_LDS(0, &mhu_a[2]);//BODYCOM0 
   v_mhu_A2    = VEC_LDS(0, &mhu_a[3]);//BODYCOM0 
   v_mhu_B1    = VEC_LDS(0, &mhu_a[9]);//BODYCOM0 
   v_mhu_B2    = VEC_LDS(0, &mhu_a[10]);//BODYCOM0 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[1]);//BODYCOM0 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[2]);//BODYCOM0 
   v_tauR_D1   = VEC_LDS(0, &tauR_a[7]);//BODYCOM0 
   v_tauR_D2   = VEC_LDS(0, &tauR_a[8]);//BODYCOM0 
//BODYCOM0 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM0 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM0 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM0 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM0 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM0 
//BODYCOM0 
//4//BODYCOM0 
/*//BODYCOM0 
   sum1[ii] = mhu_a[1]    + VM[ii]*sum1[ii];//BODYCOM0 
   sum1[ii] = mhu_a[0]    + VM[ii]*sum1[ii];//BODYCOM0 
   sum2[ii] = mhu_a[8]    + VM[ii]*sum2[ii];//BODYCOM0 
   sum3[ii] = tauRdt_a[0] + VM[ii]*sum3[ii];//BODYCOM0 
*///BODYCOM0 
//BODYCOM0 
   v_mhu_A1    = VEC_LDS(0, &mhu_a[0]);//BODYCOM0 
   v_mhu_A2    = VEC_LDS(0, &mhu_a[1]);//BODYCOM0 
   v_mhu_B1    = VEC_LDS(0, &mhu_a[8]);//BODYCOM0 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[0]);//BODYCOM0 
//BODYCOM0 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM0 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM0 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM0 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
//BODYCOM0 
/*//BODYCOM0 
   double tauRdt= sum3/sum4; //BODYCOM0 
   double mhu= sum1/sum2;//BODYCOM0 
    g[ii] +=  mhu*tauRdt - g[ii]*tauRdt;//BODYCOM0 
*///BODYCOM0 
/*//BODYCOM0 
   tauRdt[ii] = sum3[ii]/sum4[ii]; //BODYCOM0 
//BODYCOM0 
   mhu[ii]    = sum1[ii]/sum2[ii];//BODYCOM0 
//BODYCOM0 
   g[ii]     +=  mhu[ii]*tauRdt[ii] - g[ii]*tauRdt[ii];//BODYCOM0 
*///BODYCOM0 
//BODYCOM0 
// Still need to reuse the register I have above ... the code will look like gibberish because//BODYCOM0 
// the targets won't make sense.//BODYCOM0 
// A//BODYCOM0 
   v_xa      = vec_swdiv_nochk(v_sum3a, v_sum4a); //BODYCOM0 
   v_mhu_A2  = vec_swdiv_nochk(v_sum1a, v_sum2a); //BODYCOM0 
   v_mhu_B2  = vec_ld(0, &g[ii]);  //BODYCOM0 
   v_tauR_D2 = vec_sub(v_mhu_A2, v_mhu_B2);  //BODYCOM0 
   v_mhu_B2  = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);//BODYCOM0 

   
   //vec_st(v_mhu_B2, 0, &g[ii]);
   switch (misalignment) {
   case 1:
     g[ii+0] = v_mhu_B2 get [1];
     g[ii+1] = v_mhu_B2 get [2];
     g[ii+2] = v_mhu_B2 get [3];
     break;

   case 2:
     g[ii+0] = v_mhu_B2 get [2];
     g[ii+1] = v_mhu_B2 get [3];
     break;

   case 3:
     g[ii+0] = v_mhu_B2 get [3];
     break;

     default:
    //println("BAD MISALIGNMENT: %d\n", misalignment); fflush(stdout);
     assert(0);
   }

   ii = 4 - misalignment;
 }

 for (ii;ii<TTROUND(nCells,8);ii+=8)
 { 
   vdt v_xa   = vec_ld(0, &VM[ii]);  
   vdt v_xb   = vec_ld(0, &VM[ii+4]); 

//1
/*
//1
   sum1[ii] = mhu_a[7];
   sum1[ii] = mhu_a[6]    + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[14];
   sum2[ii] = mhu_a[13]   + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[6];
   sum3[ii] = tauRdt_a[5] + VM[ii]*sum3[ii];
   sum4[ii] = tauR_a[12]; 
   sum4[ii] = tauR_a[11]  + VM[ii]*sum4[ii]; 
*/

  vdt v_mhu_A1    = VEC_LDS(0, &mhu_a[6]);
  vdt v_mhu_A2    = VEC_LDS(0, &mhu_a[7]);
  vdt v_mhu_B1    = VEC_LDS(0, &mhu_a[13]);
  vdt v_mhu_B2    = VEC_LDS(0, &mhu_a[14]);
  vdt v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[5]);
  vdt v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[6]);
  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[11]);
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[12]);

  vdt v_sum1a = vec_madd(v_xa, v_mhu_A2,    v_mhu_A1);
  vdt v_sum1b = vec_madd(v_xb, v_mhu_A2,    v_mhu_A1);
  vdt v_sum2a = vec_madd(v_xa, v_mhu_B2,    v_mhu_B1);
  vdt v_sum2b = vec_madd(v_xb, v_mhu_B2,    v_mhu_B1);
  vdt v_sum3a = vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);
  vdt v_sum3b = vec_madd(v_xb, v_tauRdt_C2, v_tauRdt_C1);
  vdt v_sum4a = vec_madd(v_xa, v_tauR_D2,   v_tauR_D1);
  vdt v_sum4b = vec_madd(v_xb, v_tauR_D2,   v_tauR_D1);

//2
/*

   sum1[ii] = mhu_a[5]    + VM[ii]*sum1[ii];
   sum1[ii] = mhu_a[4]    + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[12]   + VM[ii]*sum2[ii];
   sum2[ii] = mhu_a[11]   + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[4] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[3] + VM[ii]*sum3[ii];
   sum4[ii] = tauR_a[10]  + VM[ii]*sum4[ii]; 
   sum4[ii] = tauR_a[9]   + VM[ii]*sum4[ii]; 
*/

   v_mhu_A1    = VEC_LDS(0, &mhu_a[4]);
   v_mhu_A2    = VEC_LDS(0, &mhu_a[5]);
   v_mhu_B1    = VEC_LDS(0, &mhu_a[11]);
   v_mhu_B2    = VEC_LDS(0, &mhu_a[12]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[3]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[4]);
   v_tauR_D1   = VEC_LDS(0, &tauR_a[9]);
   v_tauR_D2   = VEC_LDS(0, &tauR_a[10]);

   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);

   v_sum1b = vec_madd(v_xb, v_sum1b, v_mhu_A2);
   v_sum1b = vec_madd(v_xb, v_sum1b, v_mhu_A1);
   v_sum2b = vec_madd(v_xb, v_sum2b, v_mhu_B2);
   v_sum2b = vec_madd(v_xb, v_sum2b, v_mhu_B1); 
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   v_sum4b = vec_madd(v_xb, v_sum4b, v_tauR_D2);
   v_sum4b = vec_madd(v_xb, v_sum4b, v_tauR_D1);

//3
/*
   sum1[ii] = mhu_a[3]    + VM[ii]*sum1[ii];
   sum1[ii] = mhu_a[2]    + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[10]   + VM[ii]*sum2[ii];
   sum2[ii] = mhu_a[9]    + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[2] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[1] + VM[ii]*sum3[ii];
   sum4[ii] = tauR_a[8]   + VM[ii]*sum4[ii]; 
   sum4[ii] = tauR_a[7]   + VM[ii]*sum4[ii]; 
*/

   v_mhu_A1    = VEC_LDS(0, &mhu_a[2]);
   v_mhu_A2    = VEC_LDS(0, &mhu_a[3]);
   v_mhu_B1    = VEC_LDS(0, &mhu_a[9]);
   v_mhu_B2    = VEC_LDS(0, &mhu_a[10]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[1]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[2]);
   v_tauR_D1   = VEC_LDS(0, &tauR_a[7]);
   v_tauR_D2   = VEC_LDS(0, &tauR_a[8]);

   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);

   v_sum1b = vec_madd(v_xb, v_sum1b, v_mhu_A2);
   v_sum1b = vec_madd(v_xb, v_sum1b, v_mhu_A1);
   v_sum2b = vec_madd(v_xb, v_sum2b, v_mhu_B2);
   v_sum2b = vec_madd(v_xb, v_sum2b, v_mhu_B1); 
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   v_sum4b = vec_madd(v_xb, v_sum4b, v_tauR_D2);
   v_sum4b = vec_madd(v_xb, v_sum4b, v_tauR_D1);

//4
/*
   sum1[ii] = mhu_a[1]    + VM[ii]*sum1[ii];
   sum1[ii] = mhu_a[0]    + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[8]    + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[0] + VM[ii]*sum3[ii];
*/

   v_mhu_A1    = VEC_LDS(0, &mhu_a[0]);
   v_mhu_A2    = VEC_LDS(0, &mhu_a[1]);
   v_mhu_B1    = VEC_LDS(0, &mhu_a[8]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[0]);

   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);

   v_sum1b = vec_madd(v_xb, v_sum1b, v_mhu_A2);
   v_sum1b = vec_madd(v_xb, v_sum1b, v_mhu_A1);
   v_sum2b = vec_madd(v_xb, v_sum2b, v_mhu_B1); 
   v_sum3b = vec_madd(v_xb, v_sum3b, v_tauRdt_C1);


/*
   double tauRdt= sum3/sum4; 
   double mhu= sum1/sum2;
    g[ii] +=  mhu*tauRdt - g[ii]*tauRdt;
*/
/*
   tauRdt[ii] = sum3[ii]/sum4[ii]; 

   mhu[ii]    = sum1[ii]/sum2[ii];

   g[ii]     +=  mhu[ii]*tauRdt[ii] - g[ii]*tauRdt[ii];
*/

// Still need to reuse the register I have above ... the code will look like gibberish because
// the targets won't make sense.
// A
   v_xa      = vec_swdiv_nochk(v_sum3a, v_sum4a); 
   v_mhu_A2  = vec_swdiv_nochk(v_sum1a, v_sum2a); 
   v_mhu_B2  = vec_ld(0, &g[ii]);  
   v_tauR_D2 = vec_sub(v_mhu_A2, v_mhu_B2);  
   v_mhu_B2  = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);
   vec_st(v_mhu_B2, 0, &g[ii]);

// B
   v_xb      = vec_swdiv_nochk(v_sum3b, v_sum4b); 
   v_mhu_A1  = vec_swdiv_nochk(v_sum1b, v_sum2b); 
   v_mhu_B1  = vec_ld(0, &g[ii+4]);  
   v_tauR_D1 = vec_sub(v_mhu_A1, v_mhu_B1);  
   v_mhu_B1  = vec_madd(v_tauR_D1, v_xb, v_mhu_B1);
   vec_st(v_mhu_B1, 0, &g[ii+4]);


  } //xx
 for (ii;ii<TTROUND(nCells,4);ii+=4)
 { 
 //BODYSTART
   vdt v_xa   = vec_ld(0, &VM[ii]);  

//1
/*
//1
   sum1[ii] = mhu_a[7];
   sum1[ii] = mhu_a[6]    + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[14];
   sum2[ii] = mhu_a[13]   + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[6];
   sum3[ii] = tauRdt_a[5] + VM[ii]*sum3[ii];
   sum4[ii] = tauR_a[12]; 
   sum4[ii] = tauR_a[11]  + VM[ii]*sum4[ii]; 
*/

  vdt v_mhu_A1    = VEC_LDS(0, &mhu_a[6]);
  vdt v_mhu_A2    = VEC_LDS(0, &mhu_a[7]);
  vdt v_mhu_B1    = VEC_LDS(0, &mhu_a[13]);
  vdt v_mhu_B2    = VEC_LDS(0, &mhu_a[14]);
  vdt v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[5]);
  vdt v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[6]);
  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[11]);
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[12]);

  vdt v_sum1a = vec_madd(v_xa, v_mhu_A2,    v_mhu_A1);
  vdt v_sum2a = vec_madd(v_xa, v_mhu_B2,    v_mhu_B1);
  vdt v_sum3a = vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);
  vdt v_sum4a = vec_madd(v_xa, v_tauR_D2,   v_tauR_D1);

//2
/*

   sum1[ii] = mhu_a[5]    + VM[ii]*sum1[ii];
   sum1[ii] = mhu_a[4]    + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[12]   + VM[ii]*sum2[ii];
   sum2[ii] = mhu_a[11]   + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[4] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[3] + VM[ii]*sum3[ii];
   sum4[ii] = tauR_a[10]  + VM[ii]*sum4[ii]; 
   sum4[ii] = tauR_a[9]   + VM[ii]*sum4[ii]; 
*/

   v_mhu_A1    = VEC_LDS(0, &mhu_a[4]);
   v_mhu_A2    = VEC_LDS(0, &mhu_a[5]);
   v_mhu_B1    = VEC_LDS(0, &mhu_a[11]);
   v_mhu_B2    = VEC_LDS(0, &mhu_a[12]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[3]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[4]);
   v_tauR_D1   = VEC_LDS(0, &tauR_a[9]);
   v_tauR_D2   = VEC_LDS(0, &tauR_a[10]);

   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);

//3
/*
   sum1[ii] = mhu_a[3]    + VM[ii]*sum1[ii];
   sum1[ii] = mhu_a[2]    + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[10]   + VM[ii]*sum2[ii];
   sum2[ii] = mhu_a[9]    + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[2] + VM[ii]*sum3[ii];
   sum3[ii] = tauRdt_a[1] + VM[ii]*sum3[ii];
   sum4[ii] = tauR_a[8]   + VM[ii]*sum4[ii]; 
   sum4[ii] = tauR_a[7]   + VM[ii]*sum4[ii]; 
*/

   v_mhu_A1    = VEC_LDS(0, &mhu_a[2]);
   v_mhu_A2    = VEC_LDS(0, &mhu_a[3]);
   v_mhu_B1    = VEC_LDS(0, &mhu_a[9]);
   v_mhu_B2    = VEC_LDS(0, &mhu_a[10]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[1]);
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[2]);
   v_tauR_D1   = VEC_LDS(0, &tauR_a[7]);
   v_tauR_D2   = VEC_LDS(0, &tauR_a[8]);

   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);

//4
/*
   sum1[ii] = mhu_a[1]    + VM[ii]*sum1[ii];
   sum1[ii] = mhu_a[0]    + VM[ii]*sum1[ii];
   sum2[ii] = mhu_a[8]    + VM[ii]*sum2[ii];
   sum3[ii] = tauRdt_a[0] + VM[ii]*sum3[ii];
*/

   v_mhu_A1    = VEC_LDS(0, &mhu_a[0]);
   v_mhu_A2    = VEC_LDS(0, &mhu_a[1]);
   v_mhu_B1    = VEC_LDS(0, &mhu_a[8]);
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[0]);

   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);

/*
   double tauRdt= sum3/sum4; 
   double mhu= sum1/sum2;
    g[ii] +=  mhu*tauRdt - g[ii]*tauRdt;
*/
/*
   tauRdt[ii] = sum3[ii]/sum4[ii]; 

   mhu[ii]    = sum1[ii]/sum2[ii];

   g[ii]     +=  mhu[ii]*tauRdt[ii] - g[ii]*tauRdt[ii];
*/

// Still need to reuse the register I have above ... the code will look like gibberish because
// the targets won't make sense.
// A
   v_xa      = vec_swdiv_nochk(v_sum3a, v_sum4a); 
   v_mhu_A2  = vec_swdiv_nochk(v_sum1a, v_sum2a); 
   v_mhu_B2  = vec_ld(0, &g[ii]);  
   v_tauR_D2 = vec_sub(v_mhu_A2, v_mhu_B2);  
   v_mhu_B2  = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);
   vec_st(v_mhu_B2, 0, &g[ii]); //BODYEND
  } 

 // the elements in the end of the vector should be ignored
 int leftover = nCells - ii;
 if (leftover > 0) { 
      vdt v_xa   = vec_ld(0, &VM[ii]);  //BODYCOM1 
//BODYCOM1 
//1//BODYCOM1 
/*//BODYCOM1 
//1//BODYCOM1 
   sum1[ii] = mhu_a[7];//BODYCOM1 
   sum1[ii] = mhu_a[6]    + VM[ii]*sum1[ii];//BODYCOM1 
   sum2[ii] = mhu_a[14];//BODYCOM1 
   sum2[ii] = mhu_a[13]   + VM[ii]*sum2[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[6];//BODYCOM1 
   sum3[ii] = tauRdt_a[5] + VM[ii]*sum3[ii];//BODYCOM1 
   sum4[ii] = tauR_a[12]; //BODYCOM1 
   sum4[ii] = tauR_a[11]  + VM[ii]*sum4[ii]; //BODYCOM1 
*///BODYCOM1 
//BODYCOM1 
  vdt v_mhu_A1    = VEC_LDS(0, &mhu_a[6]);//BODYCOM1 
  vdt v_mhu_A2    = VEC_LDS(0, &mhu_a[7]);//BODYCOM1 
  vdt v_mhu_B1    = VEC_LDS(0, &mhu_a[13]);//BODYCOM1 
  vdt v_mhu_B2    = VEC_LDS(0, &mhu_a[14]);//BODYCOM1 
  vdt v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[5]);//BODYCOM1 
  vdt v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[6]);//BODYCOM1 
  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[11]);//BODYCOM1 
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[12]);//BODYCOM1 
//BODYCOM1 
  vdt v_sum1a = vec_madd(v_xa, v_mhu_A2,    v_mhu_A1);//BODYCOM1 
  vdt v_sum2a = vec_madd(v_xa, v_mhu_B2,    v_mhu_B1);//BODYCOM1 
  vdt v_sum3a = vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);//BODYCOM1 
  vdt v_sum4a = vec_madd(v_xa, v_tauR_D2,   v_tauR_D1);//BODYCOM1 
//BODYCOM1 
//2//BODYCOM1 
/*//BODYCOM1 
//BODYCOM1 
   sum1[ii] = mhu_a[5]    + VM[ii]*sum1[ii];//BODYCOM1 
   sum1[ii] = mhu_a[4]    + VM[ii]*sum1[ii];//BODYCOM1 
   sum2[ii] = mhu_a[12]   + VM[ii]*sum2[ii];//BODYCOM1 
   sum2[ii] = mhu_a[11]   + VM[ii]*sum2[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[4] + VM[ii]*sum3[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[3] + VM[ii]*sum3[ii];//BODYCOM1 
   sum4[ii] = tauR_a[10]  + VM[ii]*sum4[ii]; //BODYCOM1 
   sum4[ii] = tauR_a[9]   + VM[ii]*sum4[ii]; //BODYCOM1 
*///BODYCOM1 
//BODYCOM1 
   v_mhu_A1    = VEC_LDS(0, &mhu_a[4]);//BODYCOM1 
   v_mhu_A2    = VEC_LDS(0, &mhu_a[5]);//BODYCOM1 
   v_mhu_B1    = VEC_LDS(0, &mhu_a[11]);//BODYCOM1 
   v_mhu_B2    = VEC_LDS(0, &mhu_a[12]);//BODYCOM1 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[3]);//BODYCOM1 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[4]);//BODYCOM1 
   v_tauR_D1   = VEC_LDS(0, &tauR_a[9]);//BODYCOM1 
   v_tauR_D2   = VEC_LDS(0, &tauR_a[10]);//BODYCOM1 
//BODYCOM1 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM1 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM1 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM1 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM1 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM1 
//BODYCOM1 
//3//BODYCOM1 
/*//BODYCOM1 
   sum1[ii] = mhu_a[3]    + VM[ii]*sum1[ii];//BODYCOM1 
   sum1[ii] = mhu_a[2]    + VM[ii]*sum1[ii];//BODYCOM1 
   sum2[ii] = mhu_a[10]   + VM[ii]*sum2[ii];//BODYCOM1 
   sum2[ii] = mhu_a[9]    + VM[ii]*sum2[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[2] + VM[ii]*sum3[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[1] + VM[ii]*sum3[ii];//BODYCOM1 
   sum4[ii] = tauR_a[8]   + VM[ii]*sum4[ii]; //BODYCOM1 
   sum4[ii] = tauR_a[7]   + VM[ii]*sum4[ii]; //BODYCOM1 
*///BODYCOM1 
//BODYCOM1 
   v_mhu_A1    = VEC_LDS(0, &mhu_a[2]);//BODYCOM1 
   v_mhu_A2    = VEC_LDS(0, &mhu_a[3]);//BODYCOM1 
   v_mhu_B1    = VEC_LDS(0, &mhu_a[9]);//BODYCOM1 
   v_mhu_B2    = VEC_LDS(0, &mhu_a[10]);//BODYCOM1 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[1]);//BODYCOM1 
   v_tauRdt_C2 = VEC_LDS(0, &tauRdt_a[2]);//BODYCOM1 
   v_tauR_D1   = VEC_LDS(0, &tauR_a[7]);//BODYCOM1 
   v_tauR_D2   = VEC_LDS(0, &tauR_a[8]);//BODYCOM1 
//BODYCOM1 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM1 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM1 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM1 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM1 
   v_sum4a = vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM1 
//BODYCOM1 
//4//BODYCOM1 
/*//BODYCOM1 
   sum1[ii] = mhu_a[1]    + VM[ii]*sum1[ii];//BODYCOM1 
   sum1[ii] = mhu_a[0]    + VM[ii]*sum1[ii];//BODYCOM1 
   sum2[ii] = mhu_a[8]    + VM[ii]*sum2[ii];//BODYCOM1 
   sum3[ii] = tauRdt_a[0] + VM[ii]*sum3[ii];//BODYCOM1 
*///BODYCOM1 
//BODYCOM1 
   v_mhu_A1    = VEC_LDS(0, &mhu_a[0]);//BODYCOM1 
   v_mhu_A2    = VEC_LDS(0, &mhu_a[1]);//BODYCOM1 
   v_mhu_B1    = VEC_LDS(0, &mhu_a[8]);//BODYCOM1 
   v_tauRdt_C1 = VEC_LDS(0, &tauRdt_a[0]);//BODYCOM1 
//BODYCOM1 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM1 
   v_sum1a = vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM1 
   v_sum2a = vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM1 
   v_sum3a = vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
//BODYCOM1 
/*//BODYCOM1 
   double tauRdt= sum3/sum4; //BODYCOM1 
   double mhu= sum1/sum2;//BODYCOM1 
    g[ii] +=  mhu*tauRdt - g[ii]*tauRdt;//BODYCOM1 
*///BODYCOM1 
/*//BODYCOM1 
   tauRdt[ii] = sum3[ii]/sum4[ii]; //BODYCOM1 
//BODYCOM1 
   mhu[ii]    = sum1[ii]/sum2[ii];//BODYCOM1 
//BODYCOM1 
   g[ii]     +=  mhu[ii]*tauRdt[ii] - g[ii]*tauRdt[ii];//BODYCOM1 
*///BODYCOM1 
//BODYCOM1 
// Still need to reuse the register I have above ... the code will look like gibberish because//BODYCOM1 
// the targets won't make sense.//BODYCOM1 
// A//BODYCOM1 
   v_xa      = vec_swdiv_nochk(v_sum3a, v_sum4a); //BODYCOM1 
   v_mhu_A2  = vec_swdiv_nochk(v_sum1a, v_sum2a); //BODYCOM1 
   v_mhu_B2  = vec_ld(0, &g[ii]);  //BODYCOM1 
   v_tauR_D2 = vec_sub(v_mhu_A2, v_mhu_B2);  //BODYCOM1 
   v_mhu_B2  = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);//BODYCOM1 

   
   //printf("%d: leftover=%d\n",gateIndex,leftover); 
   //vec_st(v_mhu_B2, 0, &g[ii]);
   switch (leftover) {
   case 1:
     g[ii+0] = v_mhu_B2 get [0];
     break;

   case 2:
     g[ii+0] = v_mhu_B2 get [0];
     g[ii+1] = v_mhu_B2 get [1];
     break;

   case 3:
     g[ii+0] = v_mhu_B2 get [0];
     g[ii+1] = v_mhu_B2 get [1];
     g[ii+2] = v_mhu_B2 get [2];
     break;

     default:
    //println("BAD LEFTOVER: %d\n", leftover); fflush(stdout);
    assert(0);
   }
 }

}





void update_s1Gate_v1(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a)
{
  //update_s1Gate(dt, nCells, VM, g, mhu_a, tauR_a); return;
  typedef vector4double vdt;

int gateIndex=11;
int mhu_l=7 ;
int mhu_m= 8;
int tauR_l=11;
int tauR_m=11;
   
 double  tauRdt_a[11];


         tauRdt_a[0]  = tauR_a[0]*dt;
         tauRdt_a[1]  = tauR_a[1]*dt;
         tauRdt_a[2]  = tauR_a[2]*dt;
         tauRdt_a[3]  = tauR_a[3]*dt;
         tauRdt_a[4]  = tauR_a[4]*dt;
         tauRdt_a[5]  = tauR_a[5]*dt;
         tauRdt_a[6]  = tauR_a[6]*dt; 
         tauRdt_a[7]  = tauR_a[7]*dt;
         tauRdt_a[8]  = tauR_a[8]*dt;
         tauRdt_a[9]  = tauR_a[9]*dt;
         tauRdt_a[10] = tauR_a[10]*dt;

 int ii=0;

 // 0 aligned
 // 1: load produces (JUNK, &VM[0], &VM[1], &VM[2])
 // 2: load produces (JUNK, JUNK, &VM[0], &VM[1])
 // 3: load produces (JUNK, JUNK, JUNK, &VM[0])
 int misalignment = ((long int)&VM[0] % 32) / 8;

 // Need to allow misalignment here for now, because the workpacket is
 // subdivided by cell teype for this equation.
// assert(misalignment == 0);

 if (misalignment) { 
      vdt v_xa   = vec_ld(0, &VM[ii]);  //BODYCOM0 
//BODYCOM0 
// 1//BODYCOM0 
  vdt v_mhu_A1  = VEC_LDS(0, &mhu_a[6]);//BODYCOM0 
  vdt v_mhu_A2  = VEC_LDS(0, &mhu_a[7]);//BODYCOM0 
  vdt v_mhu_B1 = VEC_LDS(0, &mhu_a[13]);//BODYCOM0 
  vdt v_mhu_B2 = VEC_LDS(0, &mhu_a[14]);//BODYCOM0 
  vdt v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[9]);//BODYCOM0 
  vdt v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[10]);//BODYCOM0 
  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[20]);//BODYCOM0 
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[21]);//BODYCOM0 
   vdt v_sum1a =  vec_madd(v_xa, v_mhu_A2, v_mhu_A1);//BODYCOM0 
   vdt v_sum2a =  vec_madd(v_xa, v_mhu_B2, v_mhu_B1);//BODYCOM0 
   vdt v_sum3a =  vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);//BODYCOM0 
   vdt v_sum4a =  vec_madd(v_xa, v_tauR_D2, v_tauR_D1);//BODYCOM0 
// 2//BODYCOM0 
   v_mhu_A1  = VEC_LDS(0, &mhu_a[4]);//BODYCOM0 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[5]);//BODYCOM0 
   v_mhu_B1 = VEC_LDS(0, &mhu_a[11]);//BODYCOM0 
   v_mhu_B2 = VEC_LDS(0, &mhu_a[12]);//BODYCOM0 
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[7]);//BODYCOM0 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[8]);//BODYCOM0 
  v_tauR_D1   = VEC_LDS(0, &tauR_a[18]);//BODYCOM0 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[19]);//BODYCOM0 
//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM0 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM0 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM0 
// 3//BODYCOM0 
//BODYCOM0 
   v_mhu_A1  = VEC_LDS(0, &mhu_a[2]);//BODYCOM0 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[3]);//BODYCOM0 
   v_mhu_B1 = VEC_LDS(0, &mhu_a[9]);//BODYCOM0 
   v_mhu_B2 = VEC_LDS(0, &mhu_a[10]);//BODYCOM0 
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[5]);//BODYCOM0 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[6]);//BODYCOM0 
  v_tauR_D1   = VEC_LDS(0, &tauR_a[16]);//BODYCOM0 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[17]);//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM0 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM0 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM0 
// 4//BODYCOM0 
//BODYCOM0 
   v_mhu_A1  = VEC_LDS(0, &mhu_a[0]);//BODYCOM0 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[1]);//BODYCOM0 
   v_mhu_B2  = VEC_LDS(0, &mhu_a[8]);//BODYCOM0 
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[3]);//BODYCOM0 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[4]);//BODYCOM0 
  v_tauR_D1   = VEC_LDS(0, &tauR_a[14]);//BODYCOM0 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[15]);//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM0 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM0 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM0 
// 5//BODYCOM0 
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[1]);//BODYCOM0 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[2]);//BODYCOM0 
  v_tauR_D1   = VEC_LDS(0, &tauR_a[12]);//BODYCOM0 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[13]);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM0 
// 6//BODYCOM0 
//BODYCOM0 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[0]);//BODYCOM0 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[11]);//BODYCOM0 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM0 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM0 
// A//BODYCOM0 
   vdt temp1a   = vec_swdiv_nochk(v_sum3a,v_sum4a); //BODYCOM0 
   vdt temp2a   = vec_swdiv_nochk(v_sum1a,v_sum2a); //BODYCOM0 
   vdt temp3a   = vec_ld(0, &g[ii]);  //BODYCOM0 
   temp2a   = vec_sub(temp2a, temp3a);  //BODYCOM0 
   temp3a   = vec_madd(temp2a, temp1a, temp3a);//BODYCOM0 

   
   //vec_st(temp3a, 0, &g[ii]);
   switch (misalignment) {
   case 1:
     g[ii+0] = temp3a get [1];
     g[ii+1] = temp3a get [2];
     g[ii+2] = temp3a get [3];
     break;

   case 2:
     g[ii+0] = temp3a get [2];
     g[ii+1] = temp3a get [3];
     break;

   case 3:
     g[ii+0] = temp3a get [3];
     break;

     default:
    //println("BAD MISALIGNMENT: %d\n", misalignment); fflush(stdout);
     assert(0);
   }

   ii = 4 - misalignment;
 }

 for (ii;ii<TTROUND(nCells,16);ii+=16)
 { 
   vdt v_xa   = vec_ld(0, &VM[ii]);  
   vdt v_xb   = vec_ld(0, &VM[ii+4]);  
   vdt v_xc   = vec_ld(0, &VM[ii+8]);  
   vdt v_xd   = vec_ld(0, &VM[ii+12]);  
// 1
  vdt v_mhu_A1  = VEC_LDS(0, &mhu_a[6]);
  vdt v_mhu_A2  = VEC_LDS(0, &mhu_a[7]);
  vdt v_mhu_B1 = VEC_LDS(0, &mhu_a[13]);
  vdt v_mhu_B2 = VEC_LDS(0, &mhu_a[14]);
  vdt v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[9]);
  vdt v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[10]);
  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[20]);
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[21]);
   vdt v_sum1a =  vec_madd(v_xa, v_mhu_A2, v_mhu_A1);
   vdt v_sum1b =  vec_madd(v_xb, v_mhu_A2, v_mhu_A1);
   vdt v_sum1c =  vec_madd(v_xc, v_mhu_A2, v_mhu_A1);
   vdt v_sum1d =  vec_madd(v_xd, v_mhu_A2, v_mhu_A1);

   vdt v_sum2a =  vec_madd(v_xa, v_mhu_B2, v_mhu_B1);
   vdt v_sum2b =  vec_madd(v_xb, v_mhu_B2, v_mhu_B1);
   vdt v_sum2c =  vec_madd(v_xc, v_mhu_B2, v_mhu_B1);
   vdt v_sum2d =  vec_madd(v_xd, v_mhu_B2, v_mhu_B1);

   vdt v_sum3a =  vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);
   vdt v_sum3b =  vec_madd(v_xb, v_tauRdt_C2, v_tauRdt_C1);
   vdt v_sum3c =  vec_madd(v_xc, v_tauRdt_C2, v_tauRdt_C1);
   vdt v_sum3d =  vec_madd(v_xd, v_tauRdt_C2, v_tauRdt_C1);

   vdt v_sum4a =  vec_madd(v_xa, v_tauR_D2, v_tauR_D1);
   vdt v_sum4b =  vec_madd(v_xb, v_tauR_D2, v_tauR_D1);
   vdt v_sum4c =  vec_madd(v_xc, v_tauR_D2, v_tauR_D1);
   vdt v_sum4d =  vec_madd(v_xd, v_tauR_D2, v_tauR_D1);
// 2
   v_mhu_A1  = VEC_LDS(0, &mhu_a[4]);
   v_mhu_A2  = VEC_LDS(0, &mhu_a[5]);
   v_mhu_B1 = VEC_LDS(0, &mhu_a[11]);
   v_mhu_B2 = VEC_LDS(0, &mhu_a[12]);
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[7]);
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[8]);
  v_tauR_D1   = VEC_LDS(0, &tauR_a[18]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[19]);

   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A2);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A1);
   v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B2);
   v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B1); 
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D2);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D1);
   v_sum1c =  vec_madd(v_xc, v_sum1c, v_mhu_A2);
   v_sum1c =  vec_madd(v_xc, v_sum1c, v_mhu_A1);
   v_sum2c =  vec_madd(v_xc, v_sum2c, v_mhu_B2);
   v_sum2c =  vec_madd(v_xc, v_sum2c, v_mhu_B1); 
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C2);
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C1);
   v_sum4c =  vec_madd(v_xc, v_sum4c, v_tauR_D2);
   v_sum4c =  vec_madd(v_xc, v_sum4c, v_tauR_D1);
   v_sum1d =  vec_madd(v_xd, v_sum1d, v_mhu_A2);
   v_sum1d =  vec_madd(v_xd, v_sum1d, v_mhu_A1);
   v_sum2d =  vec_madd(v_xd, v_sum2d, v_mhu_B2);
   v_sum2d =  vec_madd(v_xd, v_sum2d, v_mhu_B1); 
   v_sum3d =  vec_madd(v_xd, v_sum3d, v_tauRdt_C2);
   v_sum3d =  vec_madd(v_xd, v_sum3d, v_tauRdt_C1);
   v_sum4d =  vec_madd(v_xd, v_sum4d, v_tauR_D2);
   v_sum4d =  vec_madd(v_xd, v_sum4d, v_tauR_D1);
// 3

   v_mhu_A1  = VEC_LDS(0, &mhu_a[2]);
   v_mhu_A2  = VEC_LDS(0, &mhu_a[3]);
   v_mhu_B1 = VEC_LDS(0, &mhu_a[9]);
   v_mhu_B2 = VEC_LDS(0, &mhu_a[10]);
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[5]);
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[6]);
  v_tauR_D1   = VEC_LDS(0, &tauR_a[16]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[17]);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A2);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A1);
   v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B2);
   v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B1); 
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D2);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D1);
   v_sum1c =  vec_madd(v_xc, v_sum1c, v_mhu_A2);
   v_sum1c =  vec_madd(v_xc, v_sum1c, v_mhu_A1);
   v_sum2c =  vec_madd(v_xc, v_sum2c, v_mhu_B2);
   v_sum2c =  vec_madd(v_xc, v_sum2c, v_mhu_B1); 
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C2);
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C1);
   v_sum4c =  vec_madd(v_xc, v_sum4c, v_tauR_D2);
   v_sum4c =  vec_madd(v_xc, v_sum4c, v_tauR_D1);
   v_sum1d =  vec_madd(v_xd, v_sum1d, v_mhu_A2);
   v_sum1d =  vec_madd(v_xd, v_sum1d, v_mhu_A1);
   v_sum2d =  vec_madd(v_xd, v_sum2d, v_mhu_B2);
   v_sum2d =  vec_madd(v_xd, v_sum2d, v_mhu_B1); 
   v_sum3d =  vec_madd(v_xd, v_sum3d, v_tauRdt_C2);
   v_sum3d =  vec_madd(v_xd, v_sum3d, v_tauRdt_C1);
   v_sum4d =  vec_madd(v_xd, v_sum4d, v_tauR_D2);
   v_sum4d =  vec_madd(v_xd, v_sum4d, v_tauR_D1);
// 4

   v_mhu_A1  = VEC_LDS(0, &mhu_a[0]);
   v_mhu_A2  = VEC_LDS(0, &mhu_a[1]);
   v_mhu_B2  = VEC_LDS(0, &mhu_a[8]);
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[3]);
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[4]);
  v_tauR_D1   = VEC_LDS(0, &tauR_a[14]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[15]);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A2);
   v_sum1b =  vec_madd(v_xb, v_sum1b, v_mhu_A1);
   v_sum2b =  vec_madd(v_xb, v_sum2b, v_mhu_B2);
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D2);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D1);
   v_sum1c =  vec_madd(v_xc, v_sum1c, v_mhu_A2);
   v_sum1c =  vec_madd(v_xc, v_sum1c, v_mhu_A1);
   v_sum2c =  vec_madd(v_xc, v_sum2c, v_mhu_B2);
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C2);
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C1);
   v_sum4c =  vec_madd(v_xc, v_sum4c, v_tauR_D2);
   v_sum4c =  vec_madd(v_xc, v_sum4c, v_tauR_D1);
   v_sum1d =  vec_madd(v_xd, v_sum1d, v_mhu_A2);
   v_sum1d =  vec_madd(v_xd, v_sum1d, v_mhu_A1);
   v_sum2d =  vec_madd(v_xd, v_sum2d, v_mhu_B2);
   v_sum3d =  vec_madd(v_xd, v_sum3d, v_tauRdt_C2);
   v_sum3d =  vec_madd(v_xd, v_sum3d, v_tauRdt_C1);
   v_sum4d =  vec_madd(v_xd, v_sum4d, v_tauR_D2);
   v_sum4d =  vec_madd(v_xd, v_sum4d, v_tauR_D1);
// 5

  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[1]);
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[2]);
  v_tauR_D1   = VEC_LDS(0, &tauR_a[12]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[13]);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C1);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D2);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D1);
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C2);
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C1);
   v_sum4c =  vec_madd(v_xc, v_sum4c, v_tauR_D2);
   v_sum4c =  vec_madd(v_xc, v_sum4c, v_tauR_D1);
   v_sum3d =  vec_madd(v_xd, v_sum3d, v_tauRdt_C2);
   v_sum3d =  vec_madd(v_xd, v_sum3d, v_tauRdt_C1);
   v_sum4d =  vec_madd(v_xd, v_sum4d, v_tauR_D2);
   v_sum4d =  vec_madd(v_xd, v_sum4d, v_tauR_D1);
// 6

  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[0]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[11]);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum3b =  vec_madd(v_xb, v_sum3b, v_tauRdt_C2);
   v_sum4b =  vec_madd(v_xb, v_sum4b, v_tauR_D2);
   v_sum3c =  vec_madd(v_xc, v_sum3c, v_tauRdt_C2);
   v_sum4c =  vec_madd(v_xc, v_sum4c, v_tauR_D2);
   v_sum3d =  vec_madd(v_xd, v_sum3d, v_tauRdt_C2);
   v_sum4d =  vec_madd(v_xd, v_sum4d, v_tauR_D2);
// A
   vdt temp1a   = vec_swdiv_nochk(v_sum3a,v_sum4a); 
   vdt temp2a   = vec_swdiv_nochk(v_sum1a,v_sum2a); 
   vdt temp3a   = vec_ld(0, &g[ii]);  
   temp2a   = vec_sub(temp2a, temp3a);  
   temp3a   = vec_madd(temp2a, temp1a, temp3a);
   vec_st(temp3a, 0, &g[ii]);
   /* This can also be implemented as the following, which has a higher flop rate and is more straightforward
   temp1a  = vec_swdiv_nochk(v_sum3a,v_sum4a); 
   temp2a  = vec_swdiv_nochk(v_sum1a,v_sum2a); 
   temp3a  = vec_ld(0, &g[ii]);
   temp2a   = vec_madd(temp1a, temp2a, temp3a);
   temp3a  = vec_nmadd(temp1a, temp3a, temp2a );
   vec_st(temp3a, 0, &g[ii]); */
// B
   vdt temp1b   = vec_swdiv_nochk(v_sum3b,v_sum4b); 
   vdt temp2b   = vec_swdiv_nochk(v_sum1b,v_sum2b); 
   vdt temp3b   = vec_ld(0, &g[ii+4]);  
   temp2b   = vec_sub(temp2b, temp3b);  
   temp3b   = vec_madd(temp2b, temp1b, temp3b);
   vec_st(temp3b, 0, &g[ii+4]);
// C
   vdt temp1c   = vec_swdiv_nochk(v_sum3c,v_sum4c); 
   vdt temp2c   = vec_swdiv_nochk(v_sum1c,v_sum2c); 
   vdt temp3c   = vec_ld(0, &g[ii+8]);  
   temp2c   = vec_sub(temp2c, temp3c);  
   temp3c   = vec_madd(temp2c, temp1c, temp3c);
   vec_st(temp3c, 0, &g[ii+8]);
// D
   vdt temp1d   = vec_swdiv_nochk(v_sum3d,v_sum4d); 
   vdt temp2d   = vec_swdiv_nochk(v_sum1d,v_sum2d); 
   vdt temp3d   = vec_ld(0, &g[ii+12]);  
   temp2d   = vec_sub(temp2d, temp3d);  
   temp3d   = vec_madd(temp2d, temp1d, temp3d);
   vec_st(temp3d, 0, &g[ii+12]);
 }
 
 for (ii;ii<TTROUND(nCells,4);ii+=4)
 { 
 //BODYSTART
   vdt v_xa   = vec_ld(0, &VM[ii]);  

// 1
  vdt v_mhu_A1  = VEC_LDS(0, &mhu_a[6]);
  vdt v_mhu_A2  = VEC_LDS(0, &mhu_a[7]);
  vdt v_mhu_B1 = VEC_LDS(0, &mhu_a[13]);
  vdt v_mhu_B2 = VEC_LDS(0, &mhu_a[14]);
  vdt v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[9]);
  vdt v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[10]);
  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[20]);
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[21]);
   vdt v_sum1a =  vec_madd(v_xa, v_mhu_A2, v_mhu_A1);
   vdt v_sum2a =  vec_madd(v_xa, v_mhu_B2, v_mhu_B1);
   vdt v_sum3a =  vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);
   vdt v_sum4a =  vec_madd(v_xa, v_tauR_D2, v_tauR_D1);
// 2
   v_mhu_A1  = VEC_LDS(0, &mhu_a[4]);
   v_mhu_A2  = VEC_LDS(0, &mhu_a[5]);
   v_mhu_B1 = VEC_LDS(0, &mhu_a[11]);
   v_mhu_B2 = VEC_LDS(0, &mhu_a[12]);
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[7]);
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[8]);
  v_tauR_D1   = VEC_LDS(0, &tauR_a[18]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[19]);

   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);
// 3

   v_mhu_A1  = VEC_LDS(0, &mhu_a[2]);
   v_mhu_A2  = VEC_LDS(0, &mhu_a[3]);
   v_mhu_B1 = VEC_LDS(0, &mhu_a[9]);
   v_mhu_B2 = VEC_LDS(0, &mhu_a[10]);
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[5]);
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[6]);
  v_tauR_D1   = VEC_LDS(0, &tauR_a[16]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[17]);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);
// 4

   v_mhu_A1  = VEC_LDS(0, &mhu_a[0]);
   v_mhu_A2  = VEC_LDS(0, &mhu_a[1]);
   v_mhu_B2  = VEC_LDS(0, &mhu_a[8]);
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[3]);
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[4]);
  v_tauR_D1   = VEC_LDS(0, &tauR_a[14]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[15]);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);
// 5
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[1]);
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[2]);
  v_tauR_D1   = VEC_LDS(0, &tauR_a[12]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[13]);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);
// 6

  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[0]);
  v_tauR_D2   = VEC_LDS(0, &tauR_a[11]);
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);
// A
   vdt temp1a   = vec_swdiv_nochk(v_sum3a,v_sum4a); 
   vdt temp2a   = vec_swdiv_nochk(v_sum1a,v_sum2a); 
   vdt temp3a   = vec_ld(0, &g[ii]);  
   temp2a   = vec_sub(temp2a, temp3a);  
   temp3a   = vec_madd(temp2a, temp1a, temp3a);
   vec_st(temp3a, 0, &g[ii]); //BODYEND
   /* This can also be implemented as the following, which has a higher flop rate and is more straightforward
   temp1a  = vec_swdiv_nochk(v_sum3a,v_sum4a); 
   temp2a  = vec_swdiv_nochk(v_sum1a,v_sum2a); 
   temp3a  = vec_ld(0, &g[ii]);
   temp2a   = vec_madd(temp1a, temp2a, temp3a);
   temp3a  = vec_nmadd(temp1a, temp3a, temp2a );
   vec_st(temp3a, 0, &g[ii]); */
 }

 // the elements in the end of the vector should be ignored
 int leftover = nCells - ii;
 if (leftover > 0)
 { 
    vdt v_xa   = vec_ld(0, &VM[ii]);  //BODYCOM1 
//BODYCOM1 
// 1//BODYCOM1 
  vdt v_mhu_A1  = VEC_LDS(0, &mhu_a[6]);//BODYCOM1 
  vdt v_mhu_A2  = VEC_LDS(0, &mhu_a[7]);//BODYCOM1 
  vdt v_mhu_B1 = VEC_LDS(0, &mhu_a[13]);//BODYCOM1 
  vdt v_mhu_B2 = VEC_LDS(0, &mhu_a[14]);//BODYCOM1 
  vdt v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[9]);//BODYCOM1 
  vdt v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[10]);//BODYCOM1 
  vdt v_tauR_D1   = VEC_LDS(0, &tauR_a[20]);//BODYCOM1 
  vdt v_tauR_D2   = VEC_LDS(0, &tauR_a[21]);//BODYCOM1 
   vdt v_sum1a =  vec_madd(v_xa, v_mhu_A2, v_mhu_A1);//BODYCOM1 
   vdt v_sum2a =  vec_madd(v_xa, v_mhu_B2, v_mhu_B1);//BODYCOM1 
   vdt v_sum3a =  vec_madd(v_xa, v_tauRdt_C2, v_tauRdt_C1);//BODYCOM1 
   vdt v_sum4a =  vec_madd(v_xa, v_tauR_D2, v_tauR_D1);//BODYCOM1 
// 2//BODYCOM1 
   v_mhu_A1  = VEC_LDS(0, &mhu_a[4]);//BODYCOM1 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[5]);//BODYCOM1 
   v_mhu_B1 = VEC_LDS(0, &mhu_a[11]);//BODYCOM1 
   v_mhu_B2 = VEC_LDS(0, &mhu_a[12]);//BODYCOM1 
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[7]);//BODYCOM1 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[8]);//BODYCOM1 
  v_tauR_D1   = VEC_LDS(0, &tauR_a[18]);//BODYCOM1 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[19]);//BODYCOM1 
//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM1 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM1 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM1 
// 3//BODYCOM1 
//BODYCOM1 
   v_mhu_A1  = VEC_LDS(0, &mhu_a[2]);//BODYCOM1 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[3]);//BODYCOM1 
   v_mhu_B1 = VEC_LDS(0, &mhu_a[9]);//BODYCOM1 
   v_mhu_B2 = VEC_LDS(0, &mhu_a[10]);//BODYCOM1 
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[5]);//BODYCOM1 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[6]);//BODYCOM1 
  v_tauR_D1   = VEC_LDS(0, &tauR_a[16]);//BODYCOM1 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[17]);//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM1 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM1 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B1); //BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM1 
// 4//BODYCOM1 
//BODYCOM1 
   v_mhu_A1  = VEC_LDS(0, &mhu_a[0]);//BODYCOM1 
   v_mhu_A2  = VEC_LDS(0, &mhu_a[1]);//BODYCOM1 
   v_mhu_B2  = VEC_LDS(0, &mhu_a[8]);//BODYCOM1 
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[3]);//BODYCOM1 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[4]);//BODYCOM1 
  v_tauR_D1   = VEC_LDS(0, &tauR_a[14]);//BODYCOM1 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[15]);//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A2);//BODYCOM1 
   v_sum1a =  vec_madd(v_xa, v_sum1a, v_mhu_A1);//BODYCOM1 
   v_sum2a =  vec_madd(v_xa, v_sum2a, v_mhu_B2);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM1 
// 5//BODYCOM1 
  v_tauRdt_C1   = VEC_LDS(0, &tauRdt_a[1]);//BODYCOM1 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[2]);//BODYCOM1 
  v_tauR_D1   = VEC_LDS(0, &tauR_a[12]);//BODYCOM1 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[13]);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C1);//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D1);//BODYCOM1 
// 6//BODYCOM1 
//BODYCOM1 
  v_tauRdt_C2   = VEC_LDS(0, &tauRdt_a[0]);//BODYCOM1 
  v_tauR_D2   = VEC_LDS(0, &tauR_a[11]);//BODYCOM1 
   v_sum3a =  vec_madd(v_xa, v_sum3a, v_tauRdt_C2);//BODYCOM1 
   v_sum4a =  vec_madd(v_xa, v_sum4a, v_tauR_D2);//BODYCOM1 
// A//BODYCOM1 
   vdt temp1a   = vec_swdiv_nochk(v_sum3a,v_sum4a); //BODYCOM1 
   vdt temp2a   = vec_swdiv_nochk(v_sum1a,v_sum2a); //BODYCOM1 
   vdt temp3a   = vec_ld(0, &g[ii]);  //BODYCOM1 
   temp2a   = vec_sub(temp2a, temp3a);  //BODYCOM1 
   temp3a   = vec_madd(temp2a, temp1a, temp3a);//BODYCOM1 

//   printf("%d: leftover=%d\n",gateIndex,leftover); 
   
   //vec_st(temp3a, 0, &g[ii]);
   switch (leftover) {
   case 1:
     g[ii+0] = temp3a get [0];
     break;

   case 2:
     g[ii+0] = temp3a get [0];
     g[ii+1] = temp3a get [1];
     break;

   case 3:
     g[ii+0] = temp3a get [0];
     g[ii+1] = temp3a get [1];
     g[ii+2] = temp3a get [2];
     break;

     default:
      printf("BAD LEFTOVER: %d\n", leftover); fflush(stdout);
      assert(0);
   }
 }

}
