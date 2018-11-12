#include <stdio.h>
#include "portableSIMD.h" 

typedef vector4double vdt;
#define   tauRdt_a tauR_a
#define VEC_LDS vec_lds
#define TTROUND(nc,n) ((nc)-((n)-1))

static  double exp_a[32]__attribute__((aligned(32)));
void initExp()
{
   exp_a[0]=1.0; 
   for (int i=1;i<32;i++) exp_a[i] = exp_a[i-1]/i;
}
void update_mGate_v2(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a)
{
int gateIndex=0;
int mhu_l=10;
int mhu_m=5;
int tauR_l=1;
int tauR_m=18;
int exp_l=16;

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
  vdt v_exp_a13 = VEC_LDS(0, &exp_a[13]);
  vdt v_exp_a14 = VEC_LDS(0, &exp_a[14]); 
  vdt v_exp_a15 = VEC_LDS(0, &exp_a[15]);
  vdt v_exp_a16 = VEC_LDS(0, &exp_a[16]);


  vdt v_ONE = vec_splats(1.0);

// for (int i=1;i<=exp_l;i++) exp_a[i] = exp_a[i-1]/i;


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

       vdt v_sum4  = vec_madd(v_sum3a, v_exp_a16,    v_exp_a15);
       v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a14);
       v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a13);
       v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a12);
       v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a11);
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

   //v_mhu_B2 = vec_nmsub(v_xa, v_mhu_B2, v_mhu_B2); 
   //v_mhu_B2 = vec_madd (v_xa, v_mhu_A2, v_mhu_B2); 
   v_xa        = vec_sub(v_mhu_A2, v_mhu_B2); 
   v_mhu_B2    = vec_madd(v_tauRdt_C2, v_xa, v_mhu_B2);

   vec_st(v_mhu_B2, 0, &g[ii]);

// B

   v_mhu_A1 = vec_swdiv_nochk(v_sum1b,v_sum2b); 

   v_sum4  = vec_madd(v_sum3b, v_exp_a16,    v_exp_a15);
   v_sum4  = vec_madd(v_sum3b, v_sum4,    v_exp_a14);
   v_sum4  = vec_madd(v_sum3b, v_sum4,    v_exp_a13);
   v_sum4  = vec_madd(v_sum3b, v_sum4,    v_exp_a12);
   v_sum4  = vec_madd(v_sum3b, v_sum4,    v_exp_a11);
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
 
   //v_mhu_B1 = vec_nmsub(v_xb, v_mhu_B1, v_mhu_B1); 
   //v_mhu_B1 = vec_madd (v_xb, v_mhu_A1, v_mhu_B1); 
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

       vdt v_sum4  = vec_madd(v_sum3a, v_exp_a16,    v_exp_a15);
       v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a14);
       v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a13);
       v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a12);
       v_sum4  = vec_madd(v_sum3a, v_sum4,    v_exp_a11);
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
 
   //v_mhu_B2 = vec_nmsub(v_xa, v_mhu_B2, v_mhu_B2); 
   //v_mhu_B2 = vec_madd (v_xa, v_mhu_A2, v_mhu_B2); 
   v_xa        = vec_sub(v_mhu_A2, v_mhu_B2); 
   v_mhu_B2    = vec_madd(v_tauRdt_C2, v_xa, v_mhu_B2);

   vec_st(v_mhu_B2, 0, &g[ii]); //BODYEND

  }
}

void update_hGate_v2(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a)
{

int gateIndex=1;
int mhu_l=7 ;
int mhu_m= 9;
int tauR_l=11;
int tauR_m=11;
   

 int ii=0;

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
   v_mhu_B2 = vec_nmsub(v_xa, v_mhu_B2, v_mhu_B2); 
   v_mhu_B2 = vec_madd (v_xa, v_mhu_A2, v_mhu_B2); 
   //v_tauR_D2 = vec_sub(v_mhu_A2, v_mhu_B2);  
   //v_mhu_B2       = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);
   vec_st(v_mhu_B2, 0, &g[ii]);
// B
   v_xb = vec_swdiv_nochk(v_sum3b,v_sum4b); 
   v_mhu_A1    = vec_swdiv_nochk(v_sum1b,v_sum2b); 
   v_mhu_B1   = vec_ld(0, &g[ii+4]);  
   v_mhu_B1 = vec_nmsub(v_xb, v_mhu_B1, v_mhu_B1); 
   v_mhu_B1 = vec_madd (v_xb, v_mhu_A1, v_mhu_B1); 
   //v_tauR_D2 = vec_sub(v_mhu_A1, v_mhu_B1);  
   //v_mhu_B1       = vec_madd(v_tauR_D2, v_xb, v_mhu_B1);
   vec_st(v_mhu_B1, 0, &g[ii+4]);
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
   v_mhu_B2 = vec_nmsub(v_xa, v_mhu_B2, v_mhu_B2); 
   v_mhu_B2 = vec_madd (v_xa, v_mhu_A2, v_mhu_B2); 
   //v_tauR_D2 = vec_sub(v_mhu_A2, v_mhu_B2);  
   //v_mhu_B2       = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);
   vec_st(v_mhu_B2, 0, &g[ii]); //BODYEND
 }

}





void update_jGate_v2(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a)
{

int gateIndex=2;
int mhu_l=7 ;
int mhu_m= 9;
int tauR_l=1;
int tauR_m=13;
   
 int ii=0;
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
   temp2a = vec_nmsub(v_sum3a, temp2a, temp2a); 
   temp2a = vec_madd (v_sum3a, temp1a , temp2a); 
   //temp1a = vec_sub(temp1a, temp2a);  
   //temp2a       = vec_madd(temp1a, v_sum3a, temp2a);
   vec_st(temp2a, 0, &g[ii]);
// 2nd
   vdt temp1b    = vec_swdiv_nochk(v_sum1b,v_sum2b); 
   vdt temp2b   = vec_ld(0, &g[ii+4]);  
   temp2b = vec_nmsub(v_sum3b, temp2b, temp2b); 
   temp2b = vec_madd (v_sum3b, temp1b , temp2b); 
   //temp1b = vec_sub(temp1b, temp2b);  
   //temp2b       = vec_madd(temp1b, v_sum3b, temp2b);
   vec_st(temp2b, 0, &g[ii+4]);
// 3rd
   vdt temp1c    = vec_swdiv_nochk(v_sum1c,v_sum2c); 
   vdt temp2c   = vec_ld(0, &g[ii+8]);  
   temp2c = vec_nmsub(v_sum3c, temp2c, temp2c); 
   temp2c = vec_madd (v_sum3c, temp1c , temp2c); 
   //temp1c = vec_sub(temp1c, temp2c);  
   //temp2c       = vec_madd(temp1c, v_sum3c, temp2c);
   vec_st(temp2c, 0, &g[ii+8]);
// 4th 
   vdt temp1d    = vec_swdiv_nochk(v_sum1d,v_sum2d); 
   vdt temp2d   = vec_ld(0, &g[ii+12]);  
   temp2d = vec_nmsub(v_sum3d, temp2d, temp2d); 
   temp2d = vec_madd (v_sum3d, temp1d , temp2d); 
   //temp1d = vec_sub(temp1d, temp2d);  
   //temp2d       = vec_madd(temp1d, v_sum3d, temp2d);
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
   temp2a = vec_nmsub(v_sum3a, temp2a, temp2a); 
   temp2a = vec_madd (v_sum3a, temp1a , temp2a); 
   //temp1a = vec_sub(temp1a, temp2a);  
   //temp2a       = vec_madd(temp1a, v_sum3a, temp2a);
   vec_st(temp2a, 0, &g[ii]); //BODYEND
 }
}


void update_Xr1Gate_v2(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a)
{

int gateIndex=3;
int mhu_l= 5 ;
int mhu_m= 8;
int tauR_l= 1;
int tauR_m=13;



  vdt v_mhu_a6 = VEC_LDS(0, &mhu_a[6]);
  vdt v_mhu_a7 = VEC_LDS(0, &mhu_a[7]);
  vdt v_sum4   = VEC_LDS(0, &tauR_a[13]);

  vdt v_TWO = vec_splats(2.0);

  int ii=0;

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
   v_mhu_B2 = vec_nmsub(v_xa, v_mhu_B2, v_mhu_B2); 
   v_mhu_B2 = vec_madd (v_xa, v_mhu_A2, v_mhu_B2); 
   //vdt v_tauR_D2 = vec_sub(v_mhu_A2, v_mhu_B2);  
   //v_mhu_B2      = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);
   vec_st(v_mhu_B2, 0, &g[ii]);

// B
   v_xb          = vec_swdiv_nochk(v_sum3b, v_sum4); 
   v_mhu_A1      = vec_swdiv_nochk(v_sum1b, v_sum2b); 
   v_mhu_B1      = vec_ld(0, &g[ii+4]);  
   v_mhu_B1 = vec_nmsub(v_xb, v_mhu_B1, v_mhu_B1); 
   v_mhu_B1 = vec_madd (v_xb, v_mhu_A1, v_mhu_B1); 
   //vdt v_tauR_D1 = vec_sub(v_mhu_A1, v_mhu_B1);  
   //v_mhu_B1  = vec_madd(v_tauR_D1, v_xb, v_mhu_B1);
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
   v_mhu_B2 = vec_nmsub(v_xa, v_mhu_B2, v_mhu_B2); 
   v_mhu_B2 = vec_madd (v_xa, v_mhu_A2, v_mhu_B2); 
   //vdt v_tauR_D2 = vec_sub(v_mhu_A2, v_mhu_B2);  
   //v_mhu_B2      = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);
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

}


void update_Xr2Gate_v2(double dt, int nCells,  double *VM, double *g, double *mhu_a, double *tauR_a)
{

int gateIndex=4;
int mhu_l=1 ;
int mhu_m= 10;
int tauR_l=1;
int tauR_m=10;
   



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

     double tauRdt_a0 =  tauR_a[0]; 
     double tauRdt_a1 =  tauR_a[1]; 
     double tauRdt_a2 =  tauR_a[2]; 
     double tauRdt_a3 =  tauR_a[3]; 
     double tauRdt_a4 =  tauR_a[4]; 
     double tauRdt_a5 =  tauR_a[5]; 
     double tauRdt_a6 =  tauR_a[6]; 
     double tauRdt_a7 =  tauR_a[7]; 
     double tauRdt_a8 =  tauR_a[8]; 
     double tauRdt_a9 =  tauR_a[9]; 

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

 int ii=0;
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

   //v_g = vec_nmsub(v_sum3, v_g, v_g); 
   //v_g = vec_madd (v_sum3, v_sum1, v_g); 
   vdt temp1     = vec_madd(v_sum1,   v_sum3,   v_g);
   v_g   = vec_nmsub( v_g,   v_sum3,   temp1);  // UNDID the change ... 

   //v_g_4 = vec_nmsub(v_sum3_4, v_g_4, v_g_4); 
   //v_g_4 = vec_madd (v_sum3_4, v_sum1_4, v_g_4); 
   vdt temp1_4   = vec_madd(v_sum1_4, v_sum3_4, v_g_4);
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
   //v_g = vec_nmsub(v_sum3, v_g, v_g); 
   //v_g = vec_madd (v_sum3, v_sum1, v_g); 
   vdt temp1     = vec_madd(v_sum1,   v_sum3,   v_g);
   v_g  = vec_nmsub( v_g,   v_sum3,   temp1);
   vec_st(v_g, 0, &g[ii]); //BODYEND
 }

}


void update_XsGate_v2(double dt, int nCells,  double *VM, double *g, double *mhu_a, double *tauR_a)
{
int gateIndex= 5;
int mhu_l= 5 ;
int mhu_m= 5;
int tauR_l= 6;
int tauR_m= 9;

 int ii=0;

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
   v_mhu_B2 = vec_nmsub(v_xa, v_mhu_B2, v_mhu_B2); 
   v_mhu_B2 = vec_madd (v_xa, v_mhu_A2, v_mhu_B2); 
   //v_tauR_D2 = vec_sub(v_mhu_A2,   v_mhu_B2);  
   //v_mhu_B2  = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);
   vec_st(v_mhu_B2, 0, &g[ii]);

// B
   v_xb      = vec_swdiv_nochk(v_sum3b, v_sum4b); 
   v_mhu_A1  = vec_swdiv_nochk(v_sum1b, v_sum2b); 
   v_mhu_B1  = vec_ld(0, &g[ii+4]);  
   v_mhu_B1 = vec_nmsub(v_xb, v_mhu_B1, v_mhu_B1); 
   v_mhu_B1 = vec_madd (v_xb, v_mhu_A1, v_mhu_B1); 
   //v_tauR_D1 = vec_sub(v_mhu_A1, v_mhu_B1);  
   //v_mhu_B1  = vec_madd(v_tauR_D1, v_xb, v_mhu_B1);
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
   v_mhu_B2 = vec_nmsub(v_xa, v_mhu_B2, v_mhu_B2); 
   v_mhu_B2 = vec_madd (v_xa, v_mhu_A2, v_mhu_B2); 
   //v_tauR_D2 = vec_sub(v_mhu_A2,   v_mhu_B2);  
   //v_mhu_B2  = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);
   vec_st(v_mhu_B2, 0, &g[ii]); //BODYEND
  }


}




// 4134
void update_rGate_v2(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a)
{

int gateIndex=6;
int mhu_l=5 ;
int mhu_m= 8; // Was 5
int tauR_l= 5;
int tauR_m= 7;
   
 int ii = 0;

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
   temp3a   = vec_nmsub(temp1a, temp3a, temp3a); 
   temp3a   = vec_madd (temp1a, temp2a, temp3a); 
// temp2a   = vec_sub(temp2a, temp3a);  
// temp3a   = vec_madd(temp2a, temp1a, temp3a);
   vec_st(temp3a, 0, &g[ii]);
// B
   vdt temp1b   = vec_swdiv_nochk(v_sum3b,v_sum4b); 
   vdt temp2b   = vec_swdiv_nochk(v_sum1b,v_sum2b); 
   vdt temp3b   = vec_ld(0, &g[ii+4]);  
   temp3b   = vec_nmsub(temp1b, temp3b, temp3b); 
   temp3b   = vec_madd (temp1b, temp2b, temp3b); 
// temp2b   = vec_sub(temp2b, temp3b);  
// temp3b   = vec_madd(temp2b, temp1b, temp3b);
   vec_st(temp3b, 0, &g[ii+4]);
// C
   vdt temp1c   = vec_swdiv_nochk(v_sum3c,v_sum4c); 
   vdt temp2c   = vec_swdiv_nochk(v_sum1c,v_sum2c); 
   vdt temp3c   = vec_ld(0, &g[ii+8]);  
   temp3c   = vec_nmsub(temp1c, temp3c, temp3c); 
   temp3c   = vec_madd (temp1c, temp2c, temp3c); 
// temp2c   = vec_sub(temp2c, temp3c);  
// temp3c   = vec_madd(temp2c, temp1c, temp3c);
   vec_st(temp3c, 0, &g[ii+8]);
// D
   vdt temp1d   = vec_swdiv_nochk(v_sum3d,v_sum4d); 
   vdt temp2d   = vec_swdiv_nochk(v_sum1d,v_sum2d); 
   vdt temp3d   = vec_ld(0, &g[ii+12]);  
   temp3d   = vec_nmsub(temp1d, temp3d, temp3d); 
   temp3d   = vec_madd (temp1d, temp2d, temp3d); 
// temp2d   = vec_sub(temp2d, temp3d);  
// temp3d   = vec_madd(temp2d, temp1d, temp3d);
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
   temp3a   = vec_nmsub(temp1a, temp3a, temp3a);  
   temp3a   = vec_madd (temp1a, temp2a, temp3a); 
// temp2a   = vec_sub(temp2a, temp3a);  
// temp3a   = vec_madd(temp2a, temp1a, temp3a);
   vec_st(temp3a, 0, &g[ii]); //BODYEND
 }

}
void update_dGate_v2(double dt, int nCells,  double *VM, double *g, double *mhu_a, double *tauR_a)
{

int gateIndex= 7;
int mhu_l= 5 ;
int mhu_m= 7;
int tauR_l=7 ;
int tauR_m=10;


 vdt v_mhu_a5  = VEC_LDS(0, &mhu_a[5]);
 vdt v_mhu_a6  = VEC_LDS(0, &mhu_a[6]); 
 vdt v_mhu_a10 = VEC_LDS(0, &mhu_a[10]);
 vdt v_mhu_a11 = VEC_LDS(0, &mhu_a[11]);   

 int ii=0;

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
   v_mhu_B2 = vec_nmsub(v_xa, v_mhu_B2, v_mhu_B2); 
   v_mhu_B2 = vec_madd (v_xa, v_mhu_A2, v_mhu_B2); 
   //v_tauR_D2 = vec_sub(v_mhu_A2, v_mhu_B2);  
   //v_mhu_B2  = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);
   vec_st(v_mhu_B2, 0, &g[ii]);

// B
   v_xb      = vec_swdiv_nochk(v_sum3b, v_sum4b); 
   v_mhu_A1  = vec_swdiv_nochk(v_sum1b, v_sum2b); 
   v_mhu_B1  = vec_ld(0, &g[ii+4]);  
   v_mhu_B1 = vec_nmsub(v_xb, v_mhu_B1, v_mhu_B1); 
   v_mhu_B1 = vec_madd (v_xb, v_mhu_A1, v_mhu_B1); 
   //v_tauR_D1 = vec_sub(v_mhu_A1, v_mhu_B1);  
   //v_mhu_B1  = vec_madd(v_tauR_D1, v_xb, v_mhu_B1);
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
   v_mhu_B2 = vec_nmsub(v_xa, v_mhu_B2, v_mhu_B2); 
   v_mhu_B2 = vec_madd (v_xa, v_mhu_A2, v_mhu_B2); 
   //v_tauR_D2 = vec_sub(v_mhu_A2, v_mhu_B2);  
   //v_mhu_B2  = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);
   vec_st(v_mhu_B2, 0, &g[ii]); //BODYEND
  } 
}


void update_fGate_v2(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a)
{
int gateIndex=8;
int mhu_l=5 ;
int mhu_m= 8;
int tauR_l=13;
int tauR_m=11;
   
 int ii=0;

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
   temp3a   = vec_nmsub(temp1a, temp3a, temp3a); 
   temp3a   = vec_madd (temp1a, temp2a, temp3a); 
// temp2a   = vec_madd(temp1a, temp2a, temp3a);
// temp3a  = vec_nmadd(temp1a, temp3a, temp2a );
   vec_st(temp3a, 0, &g[ii]); */
// B
   vdt temp1b   = vec_swdiv_nochk(v_sum3b,v_sum4b); 
   vdt temp2b   = vec_swdiv_nochk(v_sum1b,v_sum2b); 
   vdt temp3b   = vec_ld(0, &g[ii+4]);  
   temp3b   = vec_nmsub(temp1b, temp3b, temp3b); 
   temp3b   = vec_madd (temp1b, temp2b, temp3b); 
// temp2b   = vec_sub(temp2b, temp3b);  
// temp3b   = vec_madd(temp2b, temp1b, temp3b);
   vec_st(temp3b, 0, &g[ii+4]);
// C
   vdt temp1c   = vec_swdiv_nochk(v_sum3c,v_sum4c); 
   vdt temp2c   = vec_swdiv_nochk(v_sum1c,v_sum2c); 
   vdt temp3c   = vec_ld(0, &g[ii+8]);  
   temp3c   = vec_nmsub(temp1c, temp3c, temp3c); 
   temp3c   = vec_madd (temp1c, temp2c, temp3c); 
// temp2c   = vec_sub(temp2c, temp3c);  
// temp3c   = vec_madd(temp2c, temp1c, temp3c);
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
   temp3a   = vec_nmsub(temp1a, temp3a, temp3a); 
   temp3a   = vec_madd (temp1a, temp2a, temp3a); 
// temp2a   = vec_sub(temp2a, temp3a);  
// temp3a   = vec_madd(temp2a, temp1a, temp3a);
   vec_st(temp3a, 0, &g[ii]); //BODYEND
 }
}



void update_f2Gate_v2(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a)
{
int gateIndex= 9;
int mhu_l= 8 ;
int mhu_m= 5;
int tauR_l=11;
int tauR_m=12;

 int ii=0;
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
   v_mhu_B2 = vec_nmsub(v_xa, v_mhu_B2, v_mhu_B2); 
   v_mhu_B2 = vec_madd (v_xa, v_mhu_A2, v_mhu_B2); 
   //v_tauR_D2 = vec_sub(v_mhu_A2, v_mhu_B2);  
   //v_mhu_B2  = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);
   vec_st(v_mhu_B2, 0, &g[ii]);

// B
   v_xb      = vec_swdiv_nochk(v_sum3b, v_sum4b); 
   v_mhu_A1  = vec_swdiv_nochk(v_sum1b, v_sum2b); 
   v_mhu_B1  = vec_ld(0, &g[ii+4]);  
   v_mhu_B1 = vec_nmsub(v_xb, v_mhu_B1, v_mhu_B1); 
   v_mhu_B1 = vec_madd (v_xb, v_mhu_A1, v_mhu_B1); 
   //v_tauR_D1 = vec_sub(v_mhu_A1, v_mhu_B1);  
   //v_mhu_B1  = vec_madd(v_tauR_D1, v_xb, v_mhu_B1);
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
   v_mhu_B2 = vec_nmsub(v_xa, v_mhu_B2, v_mhu_B2); 
   v_mhu_B2 = vec_madd (v_xa, v_mhu_A2, v_mhu_B2); 
   //v_tauR_D2 = vec_sub(v_mhu_A2, v_mhu_B2);  
   //v_mhu_B2  = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);
   vec_st(v_mhu_B2, 0, &g[ii]); //BODYEND
 } 

}

void update_jLGate_v2(double dt, int nCells,  double *VM, double *g, double *mhu_a, double *tauR_a)
{
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


  vdt temp1   = vec_splats(tauR_a[0]);


 int ii=0;

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
   temp3a   = vec_nmsub(temp1, temp3a, temp3a); 
   temp3a   = vec_madd (temp1, temp2a, temp3a); 
   //temp2a   = vec_sub(temp2a, temp3a);  
   //temp3a   = vec_madd(temp2a, temp1, temp3a);
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
   temp3b   = vec_nmsub(temp1, temp3b, temp3b); 
   temp3b   = vec_madd (temp1, temp2b, temp3b); 
   //temp2b   = vec_sub(temp2b, temp3b);  
   //temp3b   = vec_madd(temp2b, temp1, temp3b);
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
   temp3a   = vec_nmsub(temp1, temp3a, temp3a); 
   temp3a   = vec_madd (temp1, temp2a, temp3a); 
   //temp2a   = vec_sub(temp2a, temp3a);  
   //temp3a   = vec_madd(temp2a, temp1, temp3a);
   vec_st(temp3a, 0, &g[ii]); //BODYEND
   /* This can also be implemented as the following, which has a higher flop rate and is more straightforward
   // temp1  = vec_swdiv_nochk(v_sum3a,v_sum4a);  ... temp1 is as we want it
   temp2a  = vec_swdiv_nochk(v_sum1a,v_sum2a); 
   temp3a  = vec_ld(0, &g[ii]);
   temp2a   = vec_madd(temp1, temp2a, temp3a);
   temp3a  = vec_nmadd(temp1, temp3a, temp2a );
   vec_st(temp3a, 0, &g[ii]); */
 }

}
void update_s0Gate_v2(double dt, int nCells,  double *VM, double *g, double *mhu_a, double *tauR_a)
{
int gateIndex= 11;
int mhu_l= 7 ;
int mhu_m= 8;
int tauR_l=6 ;
int tauR_m=7 ;

 int ii=0;
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
   v_mhu_B2 = vec_nmsub(v_xa, v_mhu_B2, v_mhu_B2); 
   v_mhu_B2 = vec_madd (v_xa, v_mhu_A2, v_mhu_B2); 
   //v_tauR_D2 = vec_sub(v_mhu_A2, v_mhu_B2);  
   //v_mhu_B2  = vec_madd(v_tauR_D2, v_xa, v_mhu_B2); 
   vec_st(v_mhu_B2, 0, &g[ii]);

// B
   v_xb      = vec_swdiv_nochk(v_sum3b, v_sum4b); 
   v_mhu_A1  = vec_swdiv_nochk(v_sum1b, v_sum2b); 
   v_mhu_B1  = vec_ld(0, &g[ii+4]);  
   v_mhu_B1 = vec_nmsub(v_xb, v_mhu_B1, v_mhu_B1); 
   v_mhu_B1 = vec_madd (v_xb, v_mhu_A1, v_mhu_B1); 
   //v_tauR_D1 = vec_sub(v_mhu_A1, v_mhu_B1);  
   //v_mhu_B1  = vec_madd(v_tauR_D1, v_xb, v_mhu_B1);
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
   v_mhu_B2 = vec_nmsub(v_xa, v_mhu_B2, v_mhu_B2); 
   v_mhu_B2 = vec_madd (v_xa, v_mhu_A2, v_mhu_B2); 
   //v_tauR_D2 = vec_sub(v_mhu_A2, v_mhu_B2);  
   //v_mhu_B2  = vec_madd(v_tauR_D2, v_xa, v_mhu_B2);
   vec_st(v_mhu_B2, 0, &g[ii]); //BODYEND
  } 
}

void update_s1Gate_v2(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a)
{
int gateIndex=11;
int mhu_l=7 ;
int mhu_m= 8;
int tauR_l=11;
int tauR_m=11;
 int ii=0;

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
//  
//A
   vdt temp1a   = vec_swdiv_nochk(v_sum3a,v_sum4a); 
   vdt temp2a   = vec_swdiv_nochk(v_sum1a,v_sum2a); 
   vdt temp3a   = vec_ld(0, &g[ii]);  
   temp3a   = vec_nmsub(temp1a, temp3a, temp3a); 
   temp3a   = vec_madd (temp1a, temp2a, temp3a); 
// temp2a   = vec_sub(temp2a, temp3a);  
// temp3a   = vec_madd(temp2a, temp1a, temp3a);
   vec_st(temp3a, 0, &g[ii]);
// B
   vdt temp1b   = vec_swdiv_nochk(v_sum3b,v_sum4b); 
   vdt temp2b   = vec_swdiv_nochk(v_sum1b,v_sum2b); 
   vdt temp3b   = vec_ld(0, &g[ii+4]);  
   temp3b   = vec_nmsub(temp1b, temp3b, temp3b); 
   temp3b   = vec_madd (temp1b, temp2b, temp3b); 
// temp2b   = vec_sub(temp2b, temp3b);  
// temp3b   = vec_madd(temp2b, temp1b, temp3b);
   vec_st(temp3b, 0, &g[ii+4]);
// C
   vdt temp1c   = vec_swdiv_nochk(v_sum3c,v_sum4c); 
   vdt temp2c   = vec_swdiv_nochk(v_sum1c,v_sum2c); 
   vdt temp3c   = vec_ld(0, &g[ii+8]);  
   temp3c   = vec_nmsub(temp1c, temp3c, temp3c); 
   temp3c   = vec_madd (temp1c, temp2c, temp3c); 
// temp2c   = vec_sub(temp2c, temp3c);  
// temp3c   = vec_madd(temp2c, temp1c, temp3c);
   vec_st(temp3c, 0, &g[ii+8]);
// D
   vdt temp1d   = vec_swdiv_nochk(v_sum3d,v_sum4d); 
   vdt temp2d   = vec_swdiv_nochk(v_sum1d,v_sum2d); 
   vdt temp3d   = vec_ld(0, &g[ii+12]);  
   temp3d   = vec_nmsub(temp1d, temp3d, temp3d); 
   temp3d   = vec_madd (temp1d, temp2d, temp3d); 
// temp2d   = vec_sub(temp2d, temp3d);  
// temp3d   = vec_madd(temp2d, temp1d, temp3d);
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
   temp3a   = vec_nmsub(temp1a, temp3a, temp3a); 
   temp3a   = vec_madd (temp1a, temp2a, temp3a); 
   //temp2a   = vec_sub(temp2a, temp3a);  
   //temp3a   = vec_madd(temp2a, temp1a, temp3a);
   vec_st(temp3a, 0, &g[ii]); //BODYEND
 }


}


