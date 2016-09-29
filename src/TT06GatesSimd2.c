#include <assert.h>
#include <stdio.h>
#include "portableSIMD.h" 

typedef vector4double vdt;
#define   tauRdt_a tauR_a
#define VEC_LDS vec_lds
#define TTROUND(nc,n) ((nc)-((n)-1))

static  int nCell0; 
static  int nCell1; 
static  int rCell0; 
static  int rCell1; 
static double *s1_mhu_a; 
static double *s1_tauR_a; 
//initSplit sets things up for the one and only thread on a node that will process a list of cells that has both s0 and s1 cells. 
//nCell0_in is the number of s0 type cells on that thread.  
//s1Mhu and s1TauR are the values for the s1 type cells. 
//
void sGateInit_v2(double *s1Mhu, double *s1TauR, int nCell_s0, int nCell_s1)
{
   s1_mhu_a=s1Mhu; 
   s1_tauR_a=s1TauR; 

   rCell0 = nCell_s0%4; 
   nCell0 = nCell_s0-rCell0; 

   rCell1 = (4-rCell0)%4; 
   nCell1 = nCell_s1-rCell1; 
}


void update_sGate_v2(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a)
{

 int ii=0;
 int jj=0;
 
 for (jj;jj<TTROUND(nCell0,8);jj+=8)
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

   ii+=8;

  } //xx
 for (jj;jj<TTROUND(nCell0,4);jj+=4)
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
   ii+=4;
  } 
// Start s0 s1    Split 
{
int gateIndex= 11;
int mhu_l= 7 ;
int mhu_m= 8;
int tauR_l=6 ;
int tauR_m=7 ;
   
 int mhu_k  = mhu_m+mhu_l-1;  
 int tauR_k = tauR_m+tauR_l-1;  
 for (int kk=0;kk<rCell0;kk++)
 { 
   double x = VM[ii];  
   double sum1=0;
   double sum2=0; 
   double sum3=0;
   double sum4=0;
   for (int j=mhu_m-1;j>=0    ;j--)sum1   =   mhu_a[j] + x*sum1;
   for (int j=mhu_k  ;j>=mhu_m;j--)sum2   =   mhu_a[j] + x*sum2;
   for (int j=tauR_m-1;j>=0    ;j--)sum3  =  tauRdt_a[j] + x*sum3;
   for (int j=tauR_k  ;j>=tauR_m;j--)sum4 =  tauR_a[j] + x*sum4; 
   double tauRdt= sum3/sum4; 
   double mhu= sum1/sum2;
    g[ii] +=  mhu*tauRdt - g[ii]*tauRdt;
    ii++;
  }
}

mhu_a  = s1_mhu_a; 
tauR_a = s1_tauR_a; 

{
int gateIndex= 11;
int mhu_l= 7 ;
int mhu_m= 8;
int tauR_l=11;
int tauR_m=11;
   
 int mhu_k  = mhu_m+mhu_l-1;  ;
 int tauR_k = tauR_m+tauR_l-1;  
 for (int kk=0;kk<rCell1;kk++)
 { 
   double x = VM[ii];  
   double sum1=0;
   double sum2=0; 
   double sum3=0;
   double sum4=0;
   for (int j=mhu_m-1;j>=0    ;j--)sum1   =   mhu_a[j] + x*sum1;
   for (int j=mhu_k  ;j>=mhu_m;j--)sum2   =   mhu_a[j] + x*sum2;
   for (int j=tauR_m-1;j>=0    ;j--)sum3  =  tauRdt_a[j] + x*sum3;
   for (int j=tauR_k  ;j>=tauR_m;j--)sum4 =  tauR_a[j] + x*sum4; 
   double tauRdt= sum3/sum4; 
   double mhu= sum1/sum2;
    g[ii] +=  mhu*tauRdt - g[ii]*tauRdt;
    ii++;
  }
}
// Finish  s0 s1    Split 

 jj=0;
 for (jj;jj<TTROUND(nCell1,16);jj+=16)
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
   ii+=16;
 }
 
 for (jj;jj<TTROUND(nCell1,4);jj+=4)
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
   ii+=4; 
 }
}
