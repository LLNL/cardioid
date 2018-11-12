
void updateFunction(void *fit, double dt, struct CellTypeParms *cellTypeParms,int nCells, double *VM, int offset, double **state, double *dVdt)
{

  typedef vector4double vdt;

  
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
//-----------------------------------------

  vdt v_5_c3  =  vec_splats(0.5*SP[14]);
  vdt v_dt_c9 =  vec_mul(vec_splats(dt), vec_splats(SP[12]));
  vdt v_dt =  vec_splats(dt);


   double c_K1  = 1;        // done
   double c_Na  = 1 *SP[20];  //done
   double c_bNa = 1 *SP[21];   //done
   double c_CaL = 1;          //done
   double c_bCa = 1*SP[5];    //done
   double c_NaCa= 1;         //done
   double c_pCa = 1*SP[6];    // done
   double c_pK  = 1;           // done

   vdt v_P_NaK   = vec_splats(1*cellTypeParms->P_NaK);
   vdt v_g_to    = vec_splats(1*cellTypeParms->g_to); 
   vdt v_g_Ks    = vec_splats(1 *cellTypeParms->g_Ks);
   vdt v_g_Kr    = vec_splats(1 *cellTypeParms->g_Kr);
   vdt v_g_NaL   = vec_splats(1 *cellTypeParms->g_NaL);

   vdt v_ONE = vec_splats(1.0);

   vdt v_v ;  
   vdt v_I_NaCa, v_itmp5, v_I_xfer ;
   vdt v_dVR;
   vdt CSX, CSX0, CSX1, CSX2, CSX3;
   for (int ii=0;ii<nCells;ii+=4) 
  {
   __dcbt(&fvX_a[0]);
   __dcbt(&fvX_a[4]);
   __dcbt(&SP[0]);
   __dcbt(&SP[4]);
   // __dcbt(&VM[ii]); //no help ... no hurt either
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
   vdt v_fv1 =  vec_madd(v_v, CSX0, CSX1);
       v_fv1 =  vec_madd(v_v, v_fv1, CSX2);
       v_fv1 =  vec_madd(v_v, v_fv1, CSX3);
   CSX = vec_ld(0, &fvX_a[4]);
   CSX0 = vec_splat(CSX,0);
   CSX1 = vec_splat(CSX,1);
   CSX2 = vec_splat(CSX,2);
   CSX3 = vec_splat(CSX,3);
       v_fv1 =  vec_madd(v_v, v_fv1, CSX0); 
       v_fv1 =  vec_madd(v_v, v_fv1, CSX1);

   vdt v_fv2 =  vec_madd(v_v, CSX2, CSX3);
   //CSX = vec_ld(0, &fvX_a[8]); Big slowdown
       v_fv2 =  vec_madd(v_v, v_fv2, vec_splats(fvX_a[8]));
       v_fv2 =  vec_madd(v_v, v_fv2, vec_splats(fvX_a[9])); 
       v_fv2 =  vec_madd(v_v, v_fv2, vec_splats(fvX_a[10]));
   // __dcbt(&_Na_i[ii]);  ... hurts
   // __dcbt(&_Ca_i[ii]);  ... hurts
       v_fv2 =  vec_madd(v_v, v_fv2, vec_splats(fvX_a[11]));
       v_fv2 =  vec_madd(v_v, v_fv2, vec_splats(fvX_a[12]));
       v_fv2 =  vec_madd(v_v, v_fv2, vec_splats(fvX_a[13]));
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

   vdt v_sigm1 = vec_swdiv_nochk(v_x1, vec_add(v_ONE, v_x1));
   vdt v_sigm2 = vec_swdiv_nochk(v_x2, vec_add(v_ONE, v_x2));
   vdt v_sigm3 = vec_swdiv_nochk(v_x3, vec_add(v_ONE, v_x3));

   vdt v_tmp = vLog4(v_states_Ca_i);

   vdt v_dV3 = vec_sub(v_v,   vec_mul(v_5_c3, v_tmp));
   //CSX = vec_ld(0, &SP[4]); 3.51 vs. 3.43
       v_dV3 = vec_sub(v_dV3, vec_splats(SP[4])); // This one in CSX was a slowdown, but spills were reduced ... 4.29 to 4.34 (not a *big* slowdown as the previous was)

   v_states_Ca_SR = vec_ld(0, &_Ca_SR[ii]);
   v_states_Ca_ss = vec_ld(0, &_Ca_ss[ii]);  

   v_I_NaCa = vec_msub(v_CUBE, v_fv1, vec_mul(v_states_Ca_i, v_fv2));
   v_I_NaCa = vec_mul(vec_splats(c_NaCa),v_I_NaCa); 

   //CSX = vec_ld(0, &SP[8]);  3.51 vs. 3.43
   vdt v_I_sum = vec_madd(vec_splats(c_bCa),  v_dV3,        vec_mul(vec_splats(c_pCa), v_sigm1));

   v_itmp5 = vec_sub(v_states_Ca_i,         v_states_Ca_SR);
   vdt v_I_up =  vec_mul(vec_splats(SP[8]), v_sigm2);
   v_itmp5 = vec_madd(vec_splats(SP[7]), v_itmp5,      v_I_up);
   v_I_xfer = vec_mul(vec_splats(SP[9]),  vec_sub(v_states_Ca_ss, v_states_Ca_i));

//    _Ca_i[ii]   += (dt*c9)*(sigm3[ii]*(I_sum[ii]-I_NaCa[ii]+I_xfer[ii]*c15-itmp5[ii]*c16));

   v_tmp = vec_msub(v_I_xfer, vec_splats(SP[10]), vec_mul(v_itmp5, vec_splats(SP[11])));
   v_tmp = vec_sub(vec_mul(v_I_sum,vec_splats(0.5)), vec_sub(v_I_NaCa, v_tmp));
   v_tmp = vec_mul(v_tmp, v_sigm3);

   vec_st(vec_madd(v_dt_c9,v_tmp,v_states_Ca_i), 0, &_Ca_i[ii]);

   v_dVR = vec_neg(v_I_sum);
   //printf("%d dVR=%14.12f %14.12f %14.12f %14.12f\n",ii+0,v_dVR get [0],v_dVR get [1],v_dVR get [2],v_dVR get [3] ); 


   vdt v_dVK_i;
//------------------------------------------
   vdt v_fv0 =  vec_madd(v_v, vec_splats(fvX_a[14]), vec_splats(fvX_a[15]));
   CSX = vec_ld(0, &fvX_a[20]);
   CSX0 = vec_splat(CSX,0);
   CSX1 = vec_splat(CSX,1);
   CSX2 = vec_splat(CSX,2);
   CSX3 = vec_splat(CSX,3);
       v_fv0 =  vec_madd(v_v, v_fv0, vec_splats(fvX_a[16])); // fvX_a[16] CSX ... same slowdown as above and reduced spills
       v_fv0 =  vec_madd(v_v, v_fv0, vec_splats(fvX_a[17]));
       v_fv0 =  vec_madd(v_v, v_fv0, vec_splats(fvX_a[18]));
       v_fv0 =  vec_madd(v_v, v_fv0, vec_splats(fvX_a[19]));
       v_fv0 =  vec_madd(v_v, v_fv0, CSX0);
       v_fv0 =  vec_madd(v_v, v_fv0, CSX1);
       v_fv0 =  vec_madd(v_v, v_fv0, CSX2); 
       v_fv0 =  vec_madd(v_v, v_fv0, CSX3);

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
   vdt v_sigm0    = vec_swdiv_nochk(v_x0, vec_add(v_ONE, v_x0));
       v_sigm0    = vec_mul(v_P_NaK, v_sigm0); 


       vdt v_logK,v_logNa;
       v_logK = vLog8(v_stateK_i, v_states_Na_i, &v_logNa);

//  double dV0 = Vm -c3*log(stateK_i get [ii]) -c5;

   vdt v_c3 = vec_splats(SP[14]);
   vdt v_dV0 = vec_sub(v_v, vec_mul(v_c3 , v_logK));
       v_dV0 = vec_sub(v_dV0, vec_splats(SP[15]));


// double dV1 = Vm -c3*log(states[Na_i])-c4;

   vdt v_dV1 = vec_sub(v_v, vec_mul(v_c3 , v_logNa));
   //CSX = vec_ld(0, &SP[16]);  ... big slowdown 4.15 to 4.35
       v_dV1 = vec_sub(v_dV1, vec_splats(SP[16]));

// double dV2 = dV0[ii]-c3*log(c2*states[Na_i]/stateK_i[ii])+c5-c6;

   //__dcbt(&fvX_a[40]); ... hurts
   vdt v_dV2 = vec_swdiv_nochk(v_states_Na_i, v_stateK_i);
       v_dV2 = vec_mul(v_dV2,vec_splats(SP[17])); 
      
       v_dV2 = vLogSeries4(v_dV2);

       v_dV2 = vec_add(vec_nmsub(v_c3, v_dV2, v_dV0), vec_splats(SP[15])); 
       v_dV2 = vec_sub(v_dV2, vec_splats(SP[18])); 
   vdt v_fv5 = vec_swdiv_nochk(v_sum5, v_sum5d);

   // CSX = vec_ld(0, &fvX_a[40]);  slowdown 4.15 to 4.33
   CSX = vec_ld(0, &fvX_a[40]); 
   CSX0 = vec_splat(CSX,0);
   CSX1 = vec_splat(CSX,1);
   CSX2 = vec_splat(CSX,2);
   CSX3 = vec_splat(CSX,3);
   vdt v_sum6 =  vec_madd(v_dV0, vec_splats(fvX_a[37]), vec_splats(fvX_a[38]));
       v_sum6 =  vec_madd(v_dV0, v_sum6, vec_splats(fvX_a[39]));
       v_sum6 =  vec_madd(v_dV0, v_sum6, CSX0);
       v_sum6 =  vec_madd(v_dV0, v_sum6, CSX1);
       v_sum6 =  vec_madd(v_dV0, v_sum6, CSX2); 
       v_sum6 =  vec_madd(v_dV0, v_sum6, CSX3);


   CSX = vec_ld(0, &fvX_a[48]); 
   CSX0 = vec_splat(CSX,0);
   CSX1 = vec_splat(CSX,1);
   CSX2 = vec_splat(CSX,2);
   CSX3 = vec_splat(CSX,3);
   vdt v_sum6d =  vec_madd(v_dV0, vec_splats(fvX_a[44]),vec_splats(fvX_a[45]));
       v_sum6d =  vec_madd(v_dV0, v_sum6d, vec_splats(fvX_a[46]));
       v_sum6d =  vec_madd(v_dV0, v_sum6d, vec_splats(fvX_a[47]));
       v_sum6d =  vec_madd(v_dV0, v_sum6d, CSX0);
       v_sum6d =  vec_madd(v_dV0, v_sum6d, CSX1);
       v_sum6d =  vec_madd(v_dV0, v_sum6d, CSX2);
       v_sum6d =  vec_madd(v_dV0, v_sum6d, CSX3);
       v_sum6d =  vec_madd(v_dV0, v_sum6d, vec_splats(fvX_a[52]));
       v_sum6d =  vec_madd(v_dV0, v_sum6d, vec_splats(fvX_a[53]));
       v_sum6d =  vec_madd(v_dV0, v_sum6d, vec_splats(fvX_a[54]));
 
   vdt v_fv6 = vec_swdiv_nochk(v_sum6,v_sum6d); 

//     double fd =  fv5[ii]  +  fv6[ii];
       
   v_fv5 = vec_mul(vec_splats(c_pK),v_fv5); 
   v_fv6 = vec_mul(vec_splats(c_K1),v_fv6); 
   vdt I_pK = vec_mul(v_fv5,v_dV0); 
   vdt I_K1 = vec_mul(v_fv6,v_dV0); 
   vdt v_fd = vec_add(v_fv5, v_fv6);     

//     double tmp0 =  (fd[ii] +  g_to*rGate[ii]*sGate[ii]+ g_Kr*c11*Xr1Gate[ii]*Xr2Gate[ii] );

   vdt v_tmp0 = vec_mul(vec_ld(0, &rGate[ii]), vec_ld(0, &sGate[ii]));
   vdt I_to = vec_mul(v_tmp0,v_dV0); 
       I_to = vec_mul(I_to,v_g_to); 
       v_tmp0 = vec_madd(v_tmp0, v_g_to, v_fd);
   v_tmp     = vec_mul(vec_ld(0, &Xr1Gate[ii]), vec_ld(0, &Xr2Gate[ii]));
   v_tmp     = vec_mul(v_tmp,v_g_Kr);                 //JNG g_Kr; 
   vdt  I_Kr = vec_mul(v_tmp,v_dV0); 
   v_tmp0  = vec_madd(v_tmp, vec_splats(SP[19]), v_tmp0);

//  double tmp1 =  (c20*CUBE(mGate[ii])*hGate[ii]*jGate[ii]+c21);

   v_tmp = vec_ld(0, &mGate[ii]);
   vdt c_mgate = vec_mul(v_tmp,   v_tmp);
       c_mgate = vec_mul(c_mgate, v_tmp);
   vdt v_tmp1 = vec_mul(c_mgate,  vec_splats(c_Na));
       v_tmp1 = vec_mul(v_tmp1,  vec_ld(0, &hGate[ii]));
   vdt I_Na = vec_mul(v_tmp1, vec_ld(0,&jGate[ii])); 
       I_Na = vec_mul(I_Na,v_dV1); 
   vdt I_bNa = vec_mul(vec_splats(c_bNa),v_dV1); 
       v_tmp1 = vec_madd(v_tmp1, vec_ld(0, &jGate[ii]), vec_splats(c_bNa));


//  double tmp2 =  g_Ks*SQ(XsGate[ii]);
   
   vdt v_tmp2 = vec_ld(0, &XsGate[ii]);
       v_tmp2 = vec_mul(v_tmp2, v_tmp2);
   vdt v_xx = v_tmp2; 
       v_tmp2 = vec_mul(v_tmp2, v_g_Ks);
   vdt I_Ks =   vec_mul(v_tmp2,v_dV2); 
  //if (ii == 0) printf("g_Ks = %15.12f %15.12f %15.12f %15.12f\n",v_g_Ks get[0], v_xx get [0], v_dV2 get [0], I_Ks get[0]); 

// double I_NaK = sigm0[ii] * fv0[ii];   

   vdt v_I_NaK = vec_mul(v_sigm0, v_fv0);

// double itmp2 = I_NaCa[ii] - 1.5*I_NaK[ii]+tmp1[ii]*dV1[ii];
  
   vdt v_itmp2 = vec_madd(v_tmp1, v_dV1, v_I_NaCa); 
       v_itmp2 = vec_nmsub(vec_splats(1.5), v_I_NaK, v_itmp2);

//  double itmp3 = I_NaK[ii] + tmp0[ii]*dV0[ii] +tmp2[ii]*dV2[ii]; 
 
   vdt v_itmp3 = vec_madd(v_tmp0, v_dV0, v_I_NaK);
       v_itmp3 = vec_madd(v_tmp2, v_dV2, v_itmp3);


// double iNaL = g_NaL*CUBE(mGate[ii])*jLGate[ii]*dV1[ii];

   vdt v_iNaL = vec_mul(v_g_NaL, c_mgate);
       v_iNaL = vec_mul(vec_ld(0, &jLGate[ii]),     v_iNaL);
       v_iNaL = vec_mul(v_dV1,                      v_iNaL);

//  states[dVK_i]     += dt*itmp3[ii];

   v_dVK_i = vec_ld(0, &_dVK_i[ii]);
   v_dVK_i = vec_madd(v_dt,v_itmp3,v_dVK_i);


// states[Na_i]    += (dt*c9)*(iNaL[ii]*c22+itmp2[ii]+2.0*I_NaCa[ii]); 
// states[Na_i]    += (dt*c9)*(iNaL[ii]+itmp2[ii]+2.0*I_NaCa[ii]);    //JNG

   v_tmp = vec_add(v_iNaL, v_itmp2);
   v_tmp = vec_madd(vec_splats(2.0),    v_I_NaCa,         v_tmp);
   v_tmp = vec_mul(v_tmp,   vec_splats(SP[12]));

   vec_st(vec_madd(v_dt,v_tmp,v_states_Na_i), 0, &_Na_i[ii]); 

//  dVR[ii] -=  iNaL[ii] + itmp2[ii]+itmp3[ii];

   v_tmp = vec_add(v_iNaL, v_itmp2);
   v_tmp = vec_add(v_tmp,  v_itmp3);

   v_dVR = vec_sub(v_dVR, v_tmp);
   //printf("%d dVR=%14.12f %14.12f %14.12f %14.12f\n",ii+0,v_dVR get [0],v_dVR get [1],v_dVR get [2],v_dVR get [3] ); 
   //if (ii == 0) printf("   SIMD: dVR= %15.12f I= %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f\n", 
     //                       v_dVR get [0], v_iNaL get [0], v_I_NaK get[0], v_I_NaCa get[0], I_Na get[0], I_bNa get [0], I_pK get [0], I_K1 get [0], I_to get [0], I_Kr get [0], I_Ks get[0]); 



//------------------------------------------

   //CSX = vec_ld(0, &fvX_a[60]); bad ... 3.44 to 3.77
   CSX = vec_ld(0, &fvX_a[56]); //xxx bad ... 3.44 to 3.73
   CSX0 = vec_splat(CSX,0);
   CSX1 = vec_splat(CSX,1);
   CSX2 = vec_splat(CSX,2);
   CSX3 = vec_splat(CSX,3);
   //CSX = vec_ld(0, &fvX_a[68]); bad ... 3.44 to 3.74
   vdt v_fv3 =  vec_madd(v_v, vec_splats(fvX_a[55]), CSX0);
       v_fv3 =  vec_madd(v_v, v_fv3, CSX1);
       v_fv3 =  vec_madd(v_v, v_fv3, CSX2);
       v_fv3 =  vec_madd(v_v, v_fv3, CSX3); 
   CSX = vec_ld(0, &fvX_a[64]); //xxx bad ... 3.44 to 3.73
   CSX0 = vec_splat(CSX,0);
   CSX1 = vec_splat(CSX,1);
   CSX2 = vec_splat(CSX,2);
   CSX3 = vec_splat(CSX,3);
       v_fv3 =  vec_madd(v_v, v_fv3, vec_splats(fvX_a[60]));
       v_fv3 =  vec_madd(v_v, v_fv3, vec_splats(fvX_a[61]));
       v_fv3 =  vec_madd(v_v, v_fv3, vec_splats(fvX_a[62]));
       v_fv3 =  vec_madd(v_v, v_fv3, vec_splats(fvX_a[63]));
       v_fv3 =  vec_madd(v_v, v_fv3, CSX0); 
       v_fv3 =  vec_madd(v_v, v_fv3, CSX1);

   vdt v_fv4 =  vec_madd(v_v, CSX2, CSX3);
       v_fv4 =  vec_madd(v_v, v_fv4, vec_splats(fvX_a[68])); 
       v_fv4 =  vec_madd(v_v, v_fv4, vec_splats(fvX_a[69]));
       v_fv4 =  vec_madd(v_v, v_fv4, vec_splats(fvX_a[70]));
       v_fv4 =  vec_madd(v_v, v_fv4, vec_splats(fvX_a[71]));
       v_fv4 =  vec_madd(v_v, v_fv4, vec_splats(fvX_a[72]));
       v_fv4 =  vec_madd(v_v, v_fv4, vec_splats(fvX_a[73])); 
       v_fv4 =  vec_madd(v_v, v_fv4, vec_splats(fvX_a[74]));

//      double x4 = SQ(states[Ca_SR]*c30); 


   vdt v_x4  = vec_mul(v_states_Ca_SR, vec_splats(SP[23]));
       v_x4  = vec_mul(v_x4, v_x4);

//    double x5 = SQ(states[Ca_SR]*c31+c32); 

   // CSX = vec_ld(0, &SP[24]);  3.55
   CSX = vec_ld(0, &SP[24]); 
   CSX0 = vec_splat(CSX,0);
   CSX1 = vec_splat(CSX,1);
   CSX2 = vec_splat(CSX,2);
   CSX3 = vec_splat(CSX,3);
   vdt v_x5  = vec_madd(v_states_Ca_SR, CSX0, CSX1);
       v_x5  = vec_mul(v_x5, v_x5);

//    double x6 = SQ(states[Ca_ss]*c33+c34); 

   vdt v_x6  = vec_madd(v_states_Ca_ss, CSX2, CSX3);
       v_x6  = vec_mul(v_x6, v_x6);
/*
     double sigm4 = sigm(x4[ii]); 
     double sigm5 = sigm(x5[ii]); 
     double sigm6 = sigm(x6[ii]); 

*/
   vdt v_sigm4 = vec_swdiv_nochk(v_x4, vec_add(v_ONE, v_x4));
   vdt v_sigm5 = vec_swdiv_nochk(v_x5, vec_add(v_ONE, v_x5));
   vdt v_sigm6 = vec_swdiv_nochk(v_x6, vec_add(v_ONE, v_x6));

//      double tmp8  = (c18+c19*sigm4[ii]); //Sigm4
   vdt v_tmp8 = vec_madd(vec_splats(SP[28]), v_sigm4, vec_splats(SP[29]));

//      double tmp9  = tmp8[ii]*states[Ca_ss]+c36; 
   vdt v_tmp9 = vec_madd(v_tmp8,  v_states_Ca_ss, vec_splats(SP[30]));

//    double I_CaL = dGate[ii]*fGate[ii]*f2Gate[ii]*states[fCass_gate]*(fv4-states[Ca_ss]*fv3);

   vdt v_fCass = vec_ld(0,_fCass+ii); 
   vdt v_I_CaL = vec_nmsub(v_fv3, v_states_Ca_ss, v_fv4);
       v_I_CaL = vec_mul(v_I_CaL, v_fCass);
       v_I_CaL = vec_mul(v_I_CaL, vec_ld(0, &f2Gate[ii]));
       v_I_CaL = vec_mul(v_I_CaL, vec_ld(0, &fGate[ii]));
       v_I_CaL = vec_mul(v_I_CaL, vec_ld(0, &dGate[ii]));
       v_I_CaL = vec_mul(v_I_CaL, vec_splats(c_CaL));

//    double sigmCa_ss =   SQ(states[Ca_ss])/(tmp8[ii]*c17 + SQ(states[Ca_ss]));
   vdt v_sq        = vec_mul(v_states_Ca_ss, v_states_Ca_ss);
   vdt v_sigmCa_ss = vec_swdiv_nochk(v_sq, vec_madd(v_tmp8, vec_splats(SP[31]), v_sq)); 

//     double itmp7 = sigmCa_ss[ii]*states[R_prime]*(states[Ca_SR] - states[Ca_ss]);
   vdt v_itmp7 = vec_sub(v_states_Ca_SR, v_states_Ca_ss);
       v_itmp7 = vec_mul(v_itmp7, vec_ld(0, &_R_prime[ii]));
       v_itmp7 = vec_mul(v_itmp7, v_sigmCa_ss);



//     states[Ca_ss]   += (dt*c9)*sigm6[ii]*(I_xfer[ii]+itmp7[ii]*c14+I_CaL[ii]*c13);   

   v_tmp = vec_madd(v_itmp7, vec_splats(SP[36]), v_I_xfer);
   v_tmp = vec_madd(v_I_CaL, vec_splats(SP[37]), v_tmp);
   v_tmp = vec_mul(v_tmp, v_sigm6);

   vec_st(vec_madd(v_dt_c9,v_tmp,v_states_Ca_ss), 0, &_Ca_ss[ii]);

//     states[Ca_SR]   += (dt*c9)*sigm5[ii]*(itmp5[ii]-c40*itmp7[ii]);
   v_tmp = vec_nmsub(v_itmp7, vec_splats(SP[38]), v_itmp5);
   v_tmp = vec_mul(v_tmp, v_sigm5);

   vec_st(vec_madd(v_dt_c9,v_tmp,v_states_Ca_SR), 0, &_Ca_SR[ii]); 
  
//     states[R_prime] += (dt*c9)*(c36 - tmp9[ii]*states[R_prime]);
   vdt v_R_prime = vec_ld(0, &_R_prime[ii]);
   v_tmp = vec_nmsub(v_tmp9, v_R_prime, vec_splats(SP[30]) );

   vec_st(vec_madd(v_dt_c9,v_tmp,v_R_prime), 0, &_R_prime[ii]); 

//     double t1 = 1.0/(1.0+SQ(20*states[Ca_ss])); 
//     double mhu = 0.600000*t1[ii]+0.4000000;
//     double tauR =    1.0/(80.0*t1[ii]+2.0);
//     states[fCass_gate] +=  dt*(mhu[ii] - states[fCass_gate])*tauR[ii]; 

#if fCassForm == TT06 
   vdt v_t1 = vec_mul(vec_splats(SP[32]), v_states_Ca_ss);
   v_t1 = vec_madd(v_t1, v_t1, v_ONE);
   v_t1 = vec_swdiv_nochk(v_ONE, v_t1);
   vdt v_mhu = vec_madd(vec_splats(SP[33]), v_t1, vec_splats(SP[34]));
   vdt v_tauR = vec_swdiv_nochk(v_ONE, vec_madd(v_t1, vec_splats(SP[35]), vec_splats(2.0))); 
#endif 

#if fCassForm == RICE
   vdt xCa_ss = vec_mul(vec_splats(1000.0),v_states_Ca_ss); 
   vdt xCa_i  = vec_mul(vec_splats(1000.0),v_states_Ca_i ); 
   vdt v_mhu = vec_madd(vec_swsqrt(xCa_i),xCa_i,v_ONE); 
       v_mhu = vec_madd(vec_swsqrt(xCa_ss),xCa_ss,v_mhu); 
   v_mhu = vec_swdiv(vec_splats(SP[33]),v_mhu); 
   v_mhu =   vec_add(vec_splats(SP[34]),v_mhu); 
   vdt v_tauR =   vec_swdiv(vec_splats(0.005),v_mhu); 
#endif
   
   v_tmp = vec_sub(v_mhu, v_fCass);
   v_tmp = vec_mul(v_tmp, v_tauR);
   vec_st(vec_madd(v_dt,v_tmp,v_fCass), 0, &_fCass[ii]);

//     dVR[ii] += I_CaL[ii]; 
   v_dVR = vec_add(v_dVR, v_I_CaL);

//    states[dVK_i] +=  dt*dVR[ii] ; 

   vec_st(vec_madd(v_dt,v_dVR,v_dVK_i), 0, &_dVK_i[ii]);

//    dVdt[ii]  = dVR[ii];
   vec_st(v_dVR,    0, &dVdt[ii]);

   }

}

