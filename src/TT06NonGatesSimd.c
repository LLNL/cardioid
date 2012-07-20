#include <math.h>
#include <assert.h>
#include "portableSIMD.h" 
#include "TT06Func.h" 
#include "TT06NonGates.h" 
#include "fastLog.h" 

#define sigm(x)   ((x)/(1+(x)))
#define logSeries(x)    (log(1+(x)) )
//#define logSeries(x) ((x)*(1.0+(x)*(-0.5+(x)/3.0)))

//#define fastLog(x) (log((x)))
//#define fastLog(x) (x)

void (*fv05Func)(double Vm, double *fv);
double (*fv6Func)(double dv);

double SP[40]__attribute__((aligned(32)));

void set_SP(struct nonGateCnst cnst)
{

 SP[0] = cnst.c26 ;
 SP[1] = cnst.c27 ;
 SP[2] = cnst.c28 ;
 SP[3] = cnst.c29 ;
 SP[4] = cnst.c8 ;
 SP[5] = cnst.c7 ;
 SP[6] = cnst.c24 ;
 SP[7] = cnst.c43 ;
 SP[8] = cnst.c44 ;
 SP[9] = cnst.c23 ;
 SP[10] = cnst.c15 ;
 SP[11] = cnst.c16 ;
 SP[12] = cnst.c9 ;
 SP[13] = cnst.c25 ;
 SP[14] = cnst.c3 ;
 SP[15] = cnst.c5 ;
 SP[16] = cnst.c4 ;
 SP[17] = cnst.c2 ;
 SP[18] = cnst.c6 ;
 SP[19] = cnst.c11 ;
 SP[20] = cnst.c20 ;
 SP[21] = cnst.c21 ;
 SP[22] = cnst.c22 ;
 SP[23] = cnst.c30 ;
 SP[24] = cnst.c31 ;
 SP[25] = cnst.c32 ;
 SP[26] = cnst.c33 ;
 SP[27] = cnst.c34 ;
 SP[28] = cnst.c19 ;
 SP[29] = cnst.c18 ;
 SP[30] = cnst.c36 ;
 SP[31] = cnst.c17 ;
 SP[32] = 20.0 ;
 SP[33] = 0.6 ;
 SP[34] = 0.4 ;
 SP[35] = 80.0 ;
 SP[36] = cnst.c14 ;
 SP[37] = cnst.c13 ;
 SP[38] = cnst.c40 ;
 SP[39] = cnst.c36 ;
}



double fvX_a[75]__attribute__((aligned(32))) = {
  6.02170007043617e-16, 2.60192856479526e-13, 6.04272220479134e-11, 1.46053164590980e-08, 2.41325169024597e-06, 1.57571330412697e-04,  
  -1.52370496566743e-12, -2.55280464729594e-10, 2.94210996176821e-08, 4.00772142804151e-07, -6.00926106818024e-04, 1.18653494080701e-01, -1.19531891610964e+01, 5.40504807620401e+02, 
  -1.69089871040937e-19, -1.57414794358693e-17, 2.10687633250732e-15, -3.07155144145422e-14, -6.15332407381053e-12, 1.86878131970215e-09, -2.93525992825876e-07, 2.87019818110585e-05, -2.24272406249374e-03, -1.45499489004044e+00, 
  5.64722239204646e-17, 2.40194725811093e-14, 4.32295749601499e-12, 4.31098192118690e-10, 2.61347101393462e-08, 9.83170991434406e-07, 2.16983094649991e-05, 2.19831706117241e-04, 
  2.24306976181154e-07, -2.65106715048650e-05, 2.00764738834124e-03, -6.61573242254081e-02, 1.00000000000000e+00,
  8.93564690058337e-12, -5.91201683806020e-09, 1.50747331544023e-06, -1.75996322715896e-04, 7.99339727741565e-03, -8.89837887994794e-03, 2.97942000558175e-01,
  1.95808469577737e-16, 2.07853354660018e-14, -4.46027161749137e-12, 8.84744201829959e-10, -5.36894894634127e-08, 3.38737191083109e-06, -7.48970504345503e-05, 1.59235097855721e-03, 2.17628774116084e-02, 2.17226317175310e-01, 1.00000000000000e+00, 
  -1.74506487282187e-20, -3.65847291145083e-18, 1.36854440783922e-16, 7.81666161976579e-14, 2.49915783972090e-12, -6.96639048514165e-10, -4.97022959453108e-08, 4.30241562974461e-06, 7.91075918239993e-04, 4.60580107291175e-02, 1.03965405517718e+00,
  1.89136344887477e-15, 4.55066993718199e-13, 1.18578052624143e-11, -5.22361846923760e-09, -3.79918615957248e-07, 3.41107297247441e-05, 6.31701538469039e-03, -7.81445921919197e-01, 2.55684105370893e+01};


void update_nonGate_v1(void *fit, double dt, struct CellTypeParms *cellTypeParms, int nCells, int *cellTypeVector, double *VM, int offset, double **state, double *dVdt)
{

  typedef vector4double vdt;


  assert(offset%4 == 0);
  
  double *_Na_i = state[Na_i]+offset;
  double *_Ca_i = state[Ca_i]+offset;
  double *_Ca_ss = state[Ca_ss]+offset;
  double *_Ca_SR = state[Ca_SR]+offset;
  double *_fCass = state[fCass]+offset;
  double *_dVK_i = state[dVK_i]+offset;
  double *_R_prime = state[R_prime]+offset;
#pragma disjoint(*_Na_i, *_Ca_i) 
#pragma disjoint(*_Na_i, *_Ca_ss) 
#pragma disjoint(*_Na_i, *_Ca_SR) 
#pragma disjoint(*_Na_i, *_fCass) 
#pragma disjoint(*_Na_i, *_dVK_i) 
#pragma disjoint(*_Na_i, *_R_prime) 
#pragma disjoint(*_Ca_i, *_Ca_ss) 
#pragma disjoint(*_Ca_i, *_Ca_SR) 
#pragma disjoint(*_Ca_i, *_fCass) 
#pragma disjoint(*_Ca_i, *_dVK_i) 
#pragma disjoint(*_Ca_i, *_R_prime) 
#pragma disjoint(*_Ca_ss, *_Ca_SR) 
#pragma disjoint(*_Ca_ss, *_fCass) 
#pragma disjoint(*_Ca_ss, *_dVK_i) 
#pragma disjoint(*_Ca_ss, *_R_prime) 
#pragma disjoint(*_Ca_SR, *_fCass) 
#pragma disjoint(*_Ca_SR, *_dVK_i) 
#pragma disjoint(*_Ca_SR, *_R_prime) 
#pragma disjoint(*_fCass, *_dVK_i) 
#pragma disjoint(*_fCass, *_R_prime) 
#pragma disjoint(*_dVK_i, *_R_prime) 
  double *f2Gate = state[f2_gate]+offset; 
  double *fGate = state[f_gate]+offset; 
  double *dGate = state[d_gate]+offset; 
  double *mGate = state[m_gate]+offset; 
  double *jGate = state[j_gate]+offset; 
  double *hGate = state[h_gate]+offset; 
  double *rGate = state[r_gate]+offset; 
  double *sGate = state[s_gate]+offset; 
  double *Xr1Gate = state[Xr1_gate]+offset; 
  double *Xr2Gate = state[Xr2_gate]+offset; 
  double *XsGate = state[Xs_gate]+offset; 
  double *jLGate = state[jL_gate]+offset; 
/* #pragma disjoint(*f2Gate, *fGate)
#pragma disjoint(*f2Gate, *dGate)
#pragma disjoint(*f2Gate, *mGate)
#pragma disjoint(*f2Gate, *jGate)
#pragma disjoint(*f2Gate, *hGate)
#pragma disjoint(*f2Gate, *rGate)
#pragma disjoint(*f2Gate, *sGate)
#pragma disjoint(*f2Gate, *Xr1Gate)
#pragma disjoint(*f2Gate, *Xr2Gate)
#pragma disjoint(*f2Gate, *XsGate)
#pragma disjoint(*f2Gate, *jLGate)
#pragma disjoint(*fGate, *dGate)
#pragma disjoint(*fGate, *mGate)
#pragma disjoint(*fGate, *jGate)
#pragma disjoint(*fGate, *hGate)
#pragma disjoint(*fGate, *rGate)
#pragma disjoint(*fGate, *sGate)
#pragma disjoint(*fGate, *Xr1Gate)
#pragma disjoint(*fGate, *Xr2Gate)
#pragma disjoint(*fGate, *XsGate)
#pragma disjoint(*fGate, *jLGate)
#pragma disjoint(*dGate, *mGate)
#pragma disjoint(*dGate, *jGate)
#pragma disjoint(*dGate, *hGate)
#pragma disjoint(*dGate, *rGate)
#pragma disjoint(*dGate, *sGate)
#pragma disjoint(*dGate, *Xr1Gate)
#pragma disjoint(*dGate, *Xr2Gate)
#pragma disjoint(*dGate, *XsGate)
#pragma disjoint(*dGate, *jLGate)
#pragma disjoint(*mGate, *jGate)
#pragma disjoint(*mGate, *hGate)
#pragma disjoint(*mGate, *rGate)
#pragma disjoint(*mGate, *sGate)
#pragma disjoint(*mGate, *Xr1Gate)
#pragma disjoint(*mGate, *Xr2Gate)
#pragma disjoint(*mGate, *XsGate)
#pragma disjoint(*mGate, *jLGate)
#pragma disjoint(*hGate, *jGate)
#pragma disjoint(*hGate, *rGate)
#pragma disjoint(*hGate, *sGate)
#pragma disjoint(*hGate, *Xr1Gate)
#pragma disjoint(*hGate, *Xr2Gate)
#pragma disjoint(*hGate, *XsGate)
#pragma disjoint(*hGate, *jLGate)
#pragma disjoint(*jGate, *rGate)
#pragma disjoint(*jGate, *sGate)
#pragma disjoint(*jGate, *Xr1Gate)
#pragma disjoint(*jGate, *Xr2Gate)
#pragma disjoint(*jGate, *XsGate)
#pragma disjoint(*jGate, *jLGate)
#pragma disjoint(*rGate, *sGate)
#pragma disjoint(*rGate, *Xr1Gate)
#pragma disjoint(*rGate, *Xr2Gate)
#pragma disjoint(*rGate, *XsGate)
#pragma disjoint(*rGate, *jLGate)
#pragma disjoint(*sGate, *Xr1Gate)
#pragma disjoint(*sGate, *Xr2Gate)
#pragma disjoint(*sGate, *XsGate)
#pragma disjoint(*sGate, *jLGate)
#pragma disjoint(*Xr1Gate, *Xr2Gate)
#pragma disjoint(*Xr1Gate, *XsGate)
#pragma disjoint(*Xr1Gate, *jLGate)
#pragma disjoint(*Xr2Gate, *XsGate)
#pragma disjoint(*Xr2Gate, *jLGate)
#pragma disjoint(*XsGate, *jLGate) too much information apparently ... slows code down */
//-----------------------------------------

  //vdt v_1     =  vec_splats(1.0);
  //vdt v_15    =  vec_splats(1.5);
  //vdt v_2     =  vec_splats(2.0);
  vdt v_5_c3  =  vec_splats(0.5*SP[14]);
  vdt v_dt_c9 =  vec_mul(vec_splats(dt), vec_splats(SP[12]));


// gather into local arrays
   int cellType = cellTypeVector[0];

   vdt v_P_NaK   = vec_splats(cellTypeParms[cellType].P_NaK);
   vdt v_g_to    = vec_splats(cellTypeParms[cellType].g_to); 
   vdt v_g_Ks    = vec_splats(cellTypeParms[cellType].g_Ks);
   vdt v_g_NaL   = vec_splats(cellTypeParms[cellType].g_NaL);


   //vdt v_v = vec_ld(0, &VM[ii]); 
   vdt v_v ; //= vec_ld(0, &VM[ii]); 
   vdt v_itmp0, v_itmp5, v_itmp6 ;
   vdt v_dVR;
   vdt CSX, CSX0, CSX1, CSX2, CSX3;
  // for (int ii=nCells-4;ii>=0;ii-=4)  ... this is faster, but try it later
  for (int ii=0;ii<nCells;ii+=4) 
  {
//   double fv[8]; 


   v_v = vec_ld(0, &VM[ii]);
   vdt v_states_Na_i;
   vdt v_states_Ca_SR;
   vdt v_states_Ca_ss;
//------------------------------------------
   CSX = vec_ld(0, &fvX_a[0]);
   CSX0 = vec_splat(CSX,0);
   CSX1 = vec_splat(CSX,1);
   CSX2 = vec_splat(CSX,2);
   CSX3 = vec_splat(CSX,3);
   //vdt v_sum1 =  vec_madd(v_v, vec_splats(fvX_a[0]), vec_splats(fvX_a[1]));
   //    v_sum1 =  vec_madd(v_v, v_sum1, vec_splats(fvX_a[2]));
   //    v_sum1 =  vec_madd(v_v, v_sum1, vec_splats(fvX_a[3]));
   vdt v_sum1 =  vec_madd(v_v, CSX0, CSX1);
       v_sum1 =  vec_madd(v_v, v_sum1, CSX2);
       v_sum1 =  vec_madd(v_v, v_sum1, CSX3);
   CSX = vec_ld(0, &fvX_a[4]);
   CSX0 = vec_splat(CSX,0);
   CSX1 = vec_splat(CSX,1);
   CSX2 = vec_splat(CSX,2);
   CSX3 = vec_splat(CSX,3);
       v_sum1 =  vec_madd(v_v, v_sum1, CSX0); 
       v_sum1 =  vec_madd(v_v, v_sum1, CSX1);
       //v_sum1 =  vec_madd(v_v, v_sum1, vec_splats(fvX_a[4])); 
       //v_sum1 =  vec_madd(v_v, v_sum1, vec_splats(fvX_a[5]));

   //vdt v_sum2 =  vec_madd(v_v, vec_splats(fvX_a[6]), vec_splats(fvX_a[7]));
   vdt v_sum2 =  vec_madd(v_v, CSX2, CSX3);
   //CSX = vec_ld(0, &fvX_a[8]); Big slowdown
       v_sum2 =  vec_madd(v_v, v_sum2, vec_splats(fvX_a[8]));
       v_sum2 =  vec_madd(v_v, v_sum2, vec_splats(fvX_a[9])); 
       v_sum2 =  vec_madd(v_v, v_sum2, vec_splats(fvX_a[10]));
       v_sum2 =  vec_madd(v_v, v_sum2, vec_splats(fvX_a[11]));
       v_sum2 =  vec_madd(v_v, v_sum2, vec_splats(fvX_a[12]));
       v_sum2 =  vec_madd(v_v, v_sum2, vec_splats(fvX_a[13]));
/*
     double x1 = _Ca_i[ii]*c26; 
     double x2 = SQ(_Ca_i[ii]*c27); 
     double x3 = SQ(_Ca_i[ii]*c28+c29); 

     double sigm1 = sigm(x1); 
     double sigm2 = sigm(x2); 
     double sigm3 = sigm(x3); 
*/

   v_states_Na_i = vec_ld(0, &_Na_i[ii]);
   vdt v_CUBE   = vec_mul(vec_mul(v_states_Na_i, v_states_Na_i),   v_states_Na_i);

   vdt v_states_Ca_i = vec_ld(0, &_Ca_i[ii]);
   // CSX = vec_ld(0, &SP[0]); 3.48 vs. 3.43 ... small slowdown
   vdt v_x1 = vec_mul(v_states_Ca_i,  vec_splats(SP[0]));
   vdt v_x2 = vec_mul(v_states_Ca_i,  vec_splats(SP[1]));
       v_x2 = vec_mul(v_x2,      v_x2);
   vdt v_x3 = vec_madd(v_states_Ca_i, vec_splats(SP[2]), vec_splats(SP[3])); 
       v_x3 = vec_mul(v_x3,      v_x3);

   vdt v_sigm1 = vec_swdiv_nochk(v_x1, vec_add(vec_splats(1.0), v_x1));
   vdt v_sigm2 = vec_swdiv_nochk(v_x2, vec_add(vec_splats(1.0), v_x2));
   vdt v_sigm3 = vec_swdiv_nochk(v_x3, vec_add(vec_splats(1.0), v_x3));
/*
   fv[1] = sum1[ii]; 
   fv[2] = sum2[ii]; 

     double fv1 = fv[1];
     double fv2 = fv[2];

     double dV3 = Vm- 0.5*c3*fastLog(_Ca_i[ii]) - c8;

     itmp0[ii] = (CUBE(_Na_i[ii])*fv1-_Ca_i[ii]*fv2); 
     double itmp4 = (c7*dV3+c24*sigm1[ii]); 
     itmp5[ii]=  (c43*(_Ca_i[ii] -  _Ca_SR[ii])+c44*sigm2[ii]);      
     itmp6[ii] = (c23*(_Ca_ss[ii] - _Ca_i[ii]));
*/

// logd4(v_states) is SLOWER
   vdt v_tmp = {fastLog(v_states_Ca_i get [0]), fastLog(v_states_Ca_i get [1]), fastLog(v_states_Ca_i get [2]), fastLog(v_states_Ca_i get [3])  };

   vdt v_dV3 = vec_sub(v_v,   vec_mul(v_5_c3, v_tmp));
   //CSX = vec_ld(0, &SP[4]); 3.51 vs. 3.43
       v_dV3 = vec_sub(v_dV3, vec_splats(SP[4])); // This one in CSX was a slowdown, but spills were reduced ... 4.29 to 4.34 (not a *big* slowdown as the previous was)

   v_states_Ca_SR = vec_ld(0, &_Ca_SR[ii]);
   v_states_Ca_ss = vec_ld(0, &_Ca_ss[ii]);  

   v_itmp0 = vec_msub(v_CUBE, v_sum1, vec_mul(v_states_Ca_i, v_sum2));

   //CSX = vec_ld(0, &SP[8]);  3.51 vs. 3.43
   vdt v_itmp4 = vec_madd(vec_splats(SP[5]),  v_dV3,        vec_mul(vec_splats(SP[6]), v_sigm1));

   v_itmp5 = vec_sub(v_states_Ca_i,         v_states_Ca_SR);
   v_itmp5 = vec_madd(vec_splats(SP[7]), v_itmp5,      vec_mul(vec_splats(SP[8]), v_sigm2));
   v_itmp6 = vec_mul(vec_splats(SP[9]),  vec_sub(v_states_Ca_ss, v_states_Ca_i));

/*
    _Ca_i[ii]   += (dt*c9)*(sigm3[ii]*(itmp4[ii]-itmp0[ii]+itmp6[ii]*c15-itmp5[ii]*c16));
*/
   v_tmp = vec_msub(v_itmp6, vec_splats(SP[10]), vec_mul(v_itmp5, vec_splats(SP[11])));
   v_tmp = vec_sub(v_itmp4, vec_sub(v_itmp0, v_tmp));
   v_tmp = vec_mul(vec_mul(v_tmp, v_sigm3), v_dt_c9);

   v_states_Ca_i     = vec_add(v_states_Ca_i, v_tmp);
   vec_st(v_states_Ca_i, 0, &_Ca_i[ii]);

   v_dVR = vec_neg(vec_mul(vec_splats(2.0), v_itmp4));


   vdt v_tmp_dVK_i;
//------------------------------------------
   vdt v_sum0 =  vec_madd(v_v, vec_splats(fvX_a[14]), vec_splats(fvX_a[15]));
   CSX = vec_ld(0, &fvX_a[20]);
   CSX0 = vec_splat(CSX,0);
   CSX1 = vec_splat(CSX,1);
   CSX2 = vec_splat(CSX,2);
   CSX3 = vec_splat(CSX,3);
       v_sum0 =  vec_madd(v_v, v_sum0, vec_splats(fvX_a[16])); // fvX_a[16] CSX ... same slowdown as above and reduced spills
       v_sum0 =  vec_madd(v_v, v_sum0, vec_splats(fvX_a[17]));
       v_sum0 =  vec_madd(v_v, v_sum0, vec_splats(fvX_a[18]));
       v_sum0 =  vec_madd(v_v, v_sum0, vec_splats(fvX_a[19]));
       v_sum0 =  vec_madd(v_v, v_sum0, CSX0);
       v_sum0 =  vec_madd(v_v, v_sum0, CSX1);
       v_sum0 =  vec_madd(v_v, v_sum0, CSX2); 
       v_sum0 =  vec_madd(v_v, v_sum0, CSX3);

   CSX = vec_ld(0, &fvX_a[24]); // Doing this with fvX_a[24] causes a *slight* slowdown ... nor from 4.15 to 4.17
   CSX0 = vec_splat(CSX,0);
   CSX1 = vec_splat(CSX,1);
   CSX2 = vec_splat(CSX,2);
   CSX3 = vec_splat(CSX,3);
   //CSX = vec_ld(0, &fvX_a[28]); // Doing this with fvX_a[28] causes a *bigger* slowdown ... nor from 4.15 to 4.26
   vdt v_sum5 =  vec_madd(v_v, CSX0, CSX1);
       v_sum5 =  vec_madd(v_v, v_sum5, CSX2);
       v_sum5 =  vec_madd(v_v, v_sum5, CSX3);
       v_sum5 =  vec_madd(v_v, v_sum5, vec_splats(fvX_a[28]));
       v_sum5 =  vec_madd(v_v, v_sum5, vec_splats(fvX_a[29]));
       v_sum5 =  vec_madd(v_v, v_sum5, vec_splats(fvX_a[30])); 
       v_sum5 =  vec_madd(v_v, v_sum5, vec_splats(fvX_a[31]));

   //CSX = vec_ld(0, &fvX_a[32]);  // Slowdown to 4.27
   vdt v_sum5d =  vec_madd(v_v, vec_splats(fvX_a[32]), vec_splats(fvX_a[33]));
       v_sum5d =  vec_madd(v_v, v_sum5d, vec_splats(fvX_a[34]));
       v_sum5d =  vec_madd(v_v, v_sum5d, vec_splats(fvX_a[35]));
       v_sum5d =  vec_madd(v_v, v_sum5d, vec_splats(fvX_a[36]));

 /*
     double stateK_i  = c9*(-Vm + states[dVK_i]);
     double x0 = states[Na_i]*c25; 
     double sigm0 = P_NaK*sigm(x0); 
 
  fv[0] = sum0[ii]; 
*/

   vdt v_stateK_i = vec_mul(vec_splats(SP[12]), vec_sub(vec_ld(0, &_dVK_i[ii]), v_v));
   vdt v_x0       = vec_mul(v_states_Na_i, vec_splats(SP[13]));
   vdt v_sigm0    = vec_swdiv_nochk(v_x0, vec_add(vec_splats(1.0), v_x0));
       v_sigm0    = vec_mul(v_P_NaK, v_sigm0); 

    v_tmp get [0] = fastLog(v_stateK_i get [0]);
    v_tmp get [1] = fastLog(v_stateK_i get [1]);
    v_tmp get [2] = fastLog(v_stateK_i get [2]);
    v_tmp get [3] = fastLog(v_stateK_i get [3]);


//  double dV0 = Vm -c3*log(stateK_i get [ii]) -c5;

   vdt v_c3 = vec_splats(SP[14]);
   vdt v_dV0 = vec_sub(v_v, vec_mul(v_c3 , v_tmp));
       v_dV0 = vec_sub(v_dV0, vec_splats(SP[15]));
   /* v_tmp get[0] = fastLog(_Na_i get [0+ii]);
   v_tmp get [1] = fastLog(_Na_i get [1+ii]);
   v_tmp get [2] = fastLog(_Na_i get [2+ii]);
   v_tmp get [3] = fastLog(_Na_i get [3+ii]); */

   v_tmp get [0] = fastLog(v_states_Na_i get [0]);
   v_tmp get [1] = fastLog(v_states_Na_i get [1]);
   v_tmp get [2] = fastLog(v_states_Na_i get [2]);
   v_tmp get [3] = fastLog(v_states_Na_i get [3]);

// double dV1 = Vm -c3*log(states[Na_i])-c4;

   vdt v_dV1 = vec_sub(v_v, vec_mul(v_c3 , v_tmp));
   //CSX = vec_ld(0, &SP[16]);  ... big slowdown 4.15 to 4.35
       v_dV1 = vec_sub(v_dV1, vec_splats(SP[16]));

// double dV2 = dV0[ii]-c3*logSeries(c2*states[Na_i]/stateK_i[ii])+c5-c6;

   vdt v_dV2 = vec_swdiv_nochk(v_states_Na_i, v_stateK_i);
       v_dV2 = vec_add(vec_splats(1.0), vec_mul(vec_splats(SP[17]), v_dV2));

       v_dV2 get [0] = fastLog(v_dV2 get [0]);
       v_dV2 get [1] = fastLog(v_dV2 get [1]);
       v_dV2 get [2] = fastLog(v_dV2 get [2]);
       v_dV2 get [3] = fastLog(v_dV2 get [3]);


       v_dV2 = vec_add(vec_nmsub(v_c3, v_dV2, v_dV0), vec_splats(SP[15])); 
       v_dV2 = vec_sub(v_dV2, vec_splats(SP[18])); 
/*
   fv[5] = sum5[ii]/sum5d[ii]; 

     double fv0 = fv[0];
     double fv5 = fv[5];
*/
   vdt v_fv5 = vec_swdiv_nochk(v_sum5, v_sum5d);

//    double fv6 = fv6Func(dV0[ii]);
//    fv6[ii]   = fv6Func(v_dV0[0]);
/*
   double sum6 = 0; for (int j=6; j>=0    ;j--) sum6  =  fv6_a[j] + dv*sum6;
   double sum6d =0; for (int j=17;j>=7;j--)     sum6d =  fv6_a[j] + dv*sum6d; 
   return sum6/sum6d; 
*/
   // CSX = vec_ld(0, &fvX_a[40]);  slowdown 4.15 to 4.33
   CSX = vec_ld(0, &fvX_a[44]); 
   CSX0 = vec_splat(CSX,0);
   CSX1 = vec_splat(CSX,1);
   CSX2 = vec_splat(CSX,2);
   CSX3 = vec_splat(CSX,3);
   vdt v_sum6 =  vec_madd(v_dV0, vec_splats(fvX_a[37]), vec_splats(fvX_a[38]));
       v_sum6 =  vec_madd(v_dV0, v_sum6, vec_splats(fvX_a[39]));
       v_sum6 =  vec_madd(v_dV0, v_sum6, vec_splats(fvX_a[40]));
       v_sum6 =  vec_madd(v_dV0, v_sum6, vec_splats(fvX_a[41]));
       v_sum6 =  vec_madd(v_dV0, v_sum6, vec_splats(fvX_a[42])); 
       v_sum6 =  vec_madd(v_dV0, v_sum6, vec_splats(fvX_a[43]));


   vdt v_sum6d =  vec_madd(v_dV0, CSX0, CSX1);
       v_sum6d =  vec_madd(v_dV0, v_sum6d, CSX2);
       v_sum6d =  vec_madd(v_dV0, v_sum6d, CSX3);
       v_sum6d =  vec_madd(v_dV0, v_sum6d, vec_splats(fvX_a[48]));
       v_sum6d =  vec_madd(v_dV0, v_sum6d, vec_splats(fvX_a[49]));
       v_sum6d =  vec_madd(v_dV0, v_sum6d, vec_splats(fvX_a[50]));
       v_sum6d =  vec_madd(v_dV0, v_sum6d, vec_splats(fvX_a[51]));
       v_sum6d =  vec_madd(v_dV0, v_sum6d, vec_splats(fvX_a[52]));
       v_sum6d =  vec_madd(v_dV0, v_sum6d, vec_splats(fvX_a[53]));
       v_sum6d =  vec_madd(v_dV0, v_sum6d, vec_splats(fvX_a[54]));
 
   vdt v_fv6 = vec_swdiv_nochk(v_sum6,v_sum6d); 

//     double fd =  fv5[ii]  +  fv6[ii];
   vdt v_fd = vec_add(v_fv5, v_fv6);

//     double tmp0 =  (fd[ii] +  g_to*rGate[ii]*sGate[ii]+ c11*Xr1Gate[ii]*Xr2Gate[ii] );

   vdt v_tmp0 = vec_mul(vec_ld(0, &rGate[ii]), vec_ld(0, &sGate[ii]));
       v_tmp0 = vec_madd(v_tmp0, v_g_to, v_fd);
   v_tmp   = vec_mul(vec_ld(0, &Xr1Gate[ii]), vec_ld(0, &Xr2Gate[ii]));
   v_tmp0  = vec_madd(v_tmp, vec_splats(SP[19]), v_tmp0);

//  double tmp1 =  (c20*CUBE(mGate[ii])*hGate[ii]*jGate[ii]+c21);

   v_tmp = vec_ld(0, &mGate[ii]);
   vdt c_mgate = vec_mul(v_tmp,   v_tmp);
       c_mgate = vec_mul(c_mgate, v_tmp);
   vdt v_tmp1 = vec_mul(c_mgate,  vec_splats(SP[20]));
       v_tmp1 = vec_mul(v_tmp1,  vec_ld(0, &hGate[ii]));
       v_tmp1 = vec_madd(v_tmp1, vec_ld(0, &jGate[ii]), vec_splats(SP[21]));


//  double tmp2 =  g_Ks*SQ(XsGate[ii]);
   
   vdt v_tmp2 = vec_ld(0, &XsGate[ii]);
       v_tmp2 = vec_mul(v_tmp2, v_tmp2);
       v_tmp2 = vec_mul(v_tmp2, v_g_Ks);

// double itmpA = sigm0[ii] * fv0[ii];   

   vdt v_itmpA = vec_mul(v_sigm0, v_sum0);

// double itmp2 = itmp0[ii] - 1.5*itmpA[ii]+tmp1[ii]*dV1[ii];
  
   vdt v_itmp2 = vec_madd(v_tmp1, v_dV1, v_itmp0); 
       v_itmp2 = vec_nmsub(vec_splats(1.5), v_itmpA, v_itmp2);

//  double itmp3 = itmpA[ii] + tmp0[ii]*dV0[ii] +tmp2[ii]*dV2[ii]; 
 
   vdt v_itmp3 = vec_madd(v_tmp0, v_dV0, v_itmpA);
       v_itmp3 = vec_madd(v_tmp2, v_dV2, v_itmp3);

// double iNaL = g_NaL*CUBE(mGate[ii])*jLGate[ii]*dV1[ii];

   vdt v_iNaL = vec_mul(v_g_NaL, c_mgate);
       v_iNaL = vec_mul(vec_ld(0, &jLGate[ii]),     v_iNaL);
       v_iNaL = vec_mul(v_dV1,                      v_iNaL);

//  states[dVK_i]     += dt*itmp3[ii];
   vdt v_dt = vec_splats(dt);
   v_tmp = vec_mul(v_dt, v_itmp3);

   v_tmp_dVK_i = vec_ld(0, &_dVK_i[ii]);
   v_tmp_dVK_i = vec_add(v_tmp_dVK_i, v_tmp);


// states[Na_i]    += (dt*c9)*(iNaL[ii]*c22+itmp2[ii]+2.0*itmp0[ii]); 
// states[Na_i]    += (dt*c9)*(iNaL[ii]+itmp2[ii]+2.0*itmp0[ii]);   //JNG

   //v_tmp = vec_madd(v_iNaL, vec_splats(SP[22]), v_itmp2);
   v_tmp = vec_add(v_iNaL, v_itmp2);    //JNG
   v_tmp = vec_madd(vec_splats(2.0),    v_itmp0,         v_tmp);
   v_tmp = vec_mul(v_tmp,   vec_mul(v_dt,   vec_splats(SP[12])));

   v_states_Na_i = vec_add(v_states_Na_i, v_tmp);
   vec_st(v_states_Na_i, 0, &_Na_i[ii]); 

//  dVR[ii] +=  iNaL[ii] - itmp2[ii]-itmp3[ii];
//   v_tmp = vec_sub(v_iNaL, v_itmp2);
//   v_tmp = vec_sub(v_tmp,  v_itmp3);
//   v_dVR = vec_add(v_dVR, v_tmp);

//  dVR[ii] -=  iNaL[ii] + itmp2[ii]+itmp3[ii];  //JNG

   v_tmp = vec_add(v_iNaL, v_itmp2);
   v_tmp = vec_add(v_tmp,  v_itmp3);

   v_dVR = vec_sub(v_dVR, v_tmp);



//------------------------------------------

/*

//  double sum3 = 0; for (int j=10;j>=0    ;j--)sum3  =  fv3_a[j]+VM[ii]*sum3;


//  double sum4 = 0; for (int j=8; j>=0    ;j--)sum4  =  fv4_a[j]+VM[ii]*sum4;

*/
   //CSX = vec_ld(0, &fvX_a[60]); bad ... 3.44 to 3.77
   // CSX = vec_ld(0, &fvX_a[64]); bad ... 3.44 to 3.73
   //CSX = vec_ld(0, &fvX_a[68]); bad ... 3.44 to 3.74
   vdt v_sum3 =  vec_madd(v_v, vec_splats(fvX_a[55]), vec_splats(fvX_a[56]));
       v_sum3 =  vec_madd(v_v, v_sum3, vec_splats(fvX_a[57]));
       v_sum3 =  vec_madd(v_v, v_sum3, vec_splats(fvX_a[58]));
       v_sum3 =  vec_madd(v_v, v_sum3, vec_splats(fvX_a[59])); 
       v_sum3 =  vec_madd(v_v, v_sum3, vec_splats(fvX_a[60]));
       v_sum3 =  vec_madd(v_v, v_sum3, vec_splats(fvX_a[61]));
       v_sum3 =  vec_madd(v_v, v_sum3, vec_splats(fvX_a[62]));
       v_sum3 =  vec_madd(v_v, v_sum3, vec_splats(fvX_a[63]));
       v_sum3 =  vec_madd(v_v, v_sum3, vec_splats(fvX_a[64])); 
       v_sum3 =  vec_madd(v_v, v_sum3, vec_splats(fvX_a[65]));

   vdt v_sum4 =  vec_madd(v_v, vec_splats(fvX_a[66]), vec_splats(fvX_a[67]));
       v_sum4 =  vec_madd(v_v, v_sum4, vec_splats(fvX_a[68])); 
       v_sum4 =  vec_madd(v_v, v_sum4, vec_splats(fvX_a[69]));
       v_sum4 =  vec_madd(v_v, v_sum4, vec_splats(fvX_a[70]));
       v_sum4 =  vec_madd(v_v, v_sum4, vec_splats(fvX_a[71]));
       v_sum4 =  vec_madd(v_v, v_sum4, vec_splats(fvX_a[72]));
       v_sum4 =  vec_madd(v_v, v_sum4, vec_splats(fvX_a[73])); 
       v_sum4 =  vec_madd(v_v, v_sum4, vec_splats(fvX_a[74]));

//      double x4 = SQ(states[Ca_SR]*c30); 


   vdt v_x4  = vec_mul(v_states_Ca_SR, vec_splats(SP[23]));
       v_x4  = vec_mul(v_x4, v_x4);

//    double x5 = SQ(states[Ca_SR]*c31+c32); 

   // CSX = vec_ld(0, &SP[24]);  3.55
   vdt v_x5  = vec_madd(v_states_Ca_SR, vec_splats(SP[24]), vec_splats(SP[25]));
       v_x5  = vec_mul(v_x5, v_x5);

//    double x6 = SQ(states[Ca_ss]*c33+c34); 

   vdt v_x6  = vec_madd(v_states_Ca_ss, vec_splats(SP[26]), vec_splats(SP[27]));
       v_x6  = vec_mul(v_x6, v_x6);
/*
     double sigm4 = sigm(x4[ii]); 
     double sigm5 = sigm(x5[ii]); 
     double sigm6 = sigm(x6[ii]); 

*/
   vdt v_sigm4 = vec_swdiv_nochk(v_x4, vec_add(vec_splats(1.0), v_x4));
   vdt v_sigm5 = vec_swdiv_nochk(v_x5, vec_add(vec_splats(1.0), v_x5));
   vdt v_sigm6 = vec_swdiv_nochk(v_x6, vec_add(vec_splats(1.0), v_x6));

//      double tmp8  = (c18+c19*sigm4[ii]); //Sigm4
   vdt v_tmp8 = vec_madd(vec_splats(SP[28]), v_sigm4, vec_splats(SP[29]));

//      double tmp9  = tmp8[ii]*states[Ca_ss]+c36; 
   vdt v_tmp9 = vec_madd(v_tmp8,  v_states_Ca_ss, vec_splats(SP[30]));

//    double itmp1 = dGate[ii]*fGate[ii]*f2Gate[ii]*states[fCass_gate]*(fv4-states[Ca_ss]*fv3);

   vdt v_itmp1 = vec_nmsub(v_sum3, v_states_Ca_ss, v_sum4);
       v_itmp1 = vec_mul(v_itmp1, vec_ld(0, &_fCass[ii]));
       v_itmp1 = vec_mul(v_itmp1, vec_ld(0, &f2Gate[ii]));
       v_itmp1 = vec_mul(v_itmp1, vec_ld(0, &fGate[ii]));
       v_itmp1 = vec_mul(v_itmp1, vec_ld(0, &dGate[ii]));

//    double sigmCa_ss =   SQ(states[Ca_ss])/(tmp8[ii]*c17 + SQ(states[Ca_ss]));
   vdt v_sq        = vec_mul(v_states_Ca_ss, v_states_Ca_ss);
   vdt v_sigmCa_ss = vec_swdiv_nochk(v_sq, vec_madd(v_tmp8, vec_splats(SP[31]), v_sq)); 

//     double itmp7 = sigmCa_ss[ii]*states[R_prime]*(states[Ca_SR] - states[Ca_ss]);
   vdt v_itmp7 = vec_sub(v_states_Ca_SR, v_states_Ca_ss);
       v_itmp7 = vec_mul(v_itmp7, vec_ld(0, &_R_prime[ii]));
       v_itmp7 = vec_mul(v_itmp7, v_sigmCa_ss);

//     double t1 = 1.0/(1.0+SQ(20*states[Ca_ss])); 
  vdt v_t1 = vec_mul(vec_splats(SP[32]), v_states_Ca_ss);
      v_t1 = vec_madd(v_t1, v_t1, vec_splats(1.0));
   //   v_t1 = vec_re(v_t1);  too much error
   v_t1 = vec_swdiv_nochk(vec_splats(1.0), v_t1);

//     double mhu = 0.600000*t1[ii]+0.4000000;
   vdt v_mhu = vec_madd(vec_splats(SP[33]), v_t1, vec_splats(SP[34]));

//     double tauR =    1.0/(80.0*t1[ii]+2.0);
   vdt v_tauR = vec_swdiv_nochk(vec_splats(1.0), vec_madd(v_t1, vec_splats(SP[35]), vec_splats(2.0))); 

//     states[Ca_ss]   += (dt*c9)*sigm6[ii]*(itmp6[ii]+itmp7[ii]*c14+itmp1[ii]*c13);   
   v_tmp = vec_madd(v_itmp7, vec_splats(SP[36]), v_itmp6);
   v_tmp = vec_madd(v_itmp1, vec_splats(SP[37]), v_tmp);
   v_tmp = vec_mul(v_tmp, v_sigm6);
   v_tmp = vec_mul(v_tmp, v_dt_c9);

   v_states_Ca_ss = vec_add(v_states_Ca_ss, v_tmp);
   vec_st(v_states_Ca_ss, 0, &_Ca_ss[ii]);

//     states[Ca_SR]   += (dt*c9)*sigm5[ii]*(itmp5[ii]-c40*itmp7[ii]);
   v_tmp = vec_nmsub(v_itmp7, vec_splats(SP[38]), v_itmp5);
   v_tmp = vec_mul(v_tmp, v_sigm5);
   v_tmp = vec_mul(v_tmp, v_dt_c9);

   v_states_Ca_SR = vec_add(v_states_Ca_SR, v_tmp);
   vec_st(v_states_Ca_SR, 0, &_Ca_SR[ii]); 
  
//     states[R_prime] += (dt*c9)*(c36 - tmp9[ii]*states[R_prime]);
   v_tmp = vec_nmsub(v_tmp9, vec_ld(0, &_R_prime[ii]), vec_splats(SP[30]) );
   v_tmp = vec_mul(v_tmp, v_dt_c9);

   vdt v_tmp_R_prime = vec_ld(0, &_R_prime[ii]);
       v_tmp_R_prime = vec_add(v_tmp_R_prime, v_tmp);
   vec_st(v_tmp_R_prime, 0, &_R_prime[ii]); 

//     states[fCass_gate] +=  dt*(mhu[ii] - states[fCass_gate])*tauR[ii]; 
   v_tmp = vec_sub(v_mhu, vec_ld(0, &_fCass[ii]));
   v_tmp = vec_mul(v_tmp, v_tauR);
   v_dt = vec_splats(dt);
   v_tmp = vec_mul(v_tmp, v_dt);

   vdt v_tmp_fCass_gate = vec_ld(0, &_fCass[ii]);
       v_tmp_fCass_gate = vec_add(v_tmp_fCass_gate, v_tmp);
   vec_st(v_tmp_fCass_gate, 0, &_fCass[ii]);

//     dVR[ii] += itmp1[ii]; 
   v_dVR = vec_add(v_dVR, v_itmp1);

//    states[dVK_i] +=  dt*dVR[ii] ; 

   v_tmp = vec_mul(v_dVR, vec_splats(dt));

   v_tmp_dVK_i = vec_add(v_tmp_dVK_i, v_tmp);
   vec_st(v_tmp_dVK_i, 0, &_dVK_i[ii]);

//    dVdt[ii]  = dVR[ii];
   vec_st(v_dVR,    0, &dVdt[ii]);

   }

}

