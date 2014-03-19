#include <math.h>
#include "TT06Func.h" 
#include "TT06NonGates.h" 
#define sigm(x)   ((x)/(1.0+(x)) )


static struct nonGateCnst cnst;
static  double f1,f2,f3,f4,f5,f6,f7,f7a,f9,f9a,f10,c12; 
void initNonGateCnst()
{
 double P_NaK,g_Ks[3],g_to[3]; 
 double g_NaL =0.0; 
 double c1,c2,c3,c4,c5,c6,c7,c8,c9;
 double c10,c11,c13,c14,c15,c16,c17,c18,c19;
 double c20,c21,c22,c23,c24,c25,c26,c27,c28,c29;
 double c30,c31,c32,c33,c34,c36,c40,c43,c44;
   /*
    * pcnst[0] is R in component membrane (joule_per_mole_kelvin).
    * pcnst[1] is T in component membrane (kelvin).
    * pcnst[2] is F in component membrane (coulomb_per_millimole).
    * pcnst[3] is Cm in component membrane (microF).
    * pcnst[4] is V_c in component membrane (micrometre3).
    * pcnst[5] is stim_start in component membrane (millisecond).
    * pcnst[6] is stim_period in component membrane (millisecond).
    * pcnst[7] is stim_duration in component membrane (millisecond).
    * pcnst[8] is stim_amplitude in component membrane (picoA_per_picoF).
    * pcnst[9] is P_kna in component reversal_potentials (dimensionless).
    * pcnst[10] is K_o in component potassium_dynamics (millimolar).
    * pcnst[11] is Na_o in component sodium_dynamics (millimolar).
    * pcnst[12] is Ca_o in component calcium_dynamics (millimolar).
    * pcnst[13] is g_K1 in component inward_rectifier_potassium_current (nanoS_per_picoF).
    * pcnst[14] is g_Kr in component rapid_time_dependent_potassium_current (nanoS_per_picoF).
    * pcnst[15] is g_Ks in component slow_time_dependent_potassium_current (nanoS_per_picoF).
    * pcnst[16] is g_Na in component fast_sodium_current (nanoS_per_picoF).
    * pcnst[17] is g_bna in component sodium_background_current (nanoS_per_picoF).
    * pcnst[18] is g_CaL in component L_type_Ca_current (litre_per_farad_second).
    * pcnst[19] is g_bca in component calcium_background_current (nanoS_per_picoF).
    * pcnst[20] is g_to in component transient_outward_current (nanoS_per_picoF).
    * pcnst[21] is P_NaK in component sodium_potassium_pump_current (picoA_per_picoF).
    * pcnst[22] is K_mk in component sodium_potassium_pump_current (millimolar).
    * pcnst[23] is K_mNa in component sodium_potassium_pump_current (millimolar).
    * pcnst[24] is K_NaCa in component sodium_calcium_exchanger_current (picoA_per_picoF).
    * pcnst[25] is K_sat in component sodium_calcium_exchanger_current (dimensionless).
    * pcnst[26] is alpha in component sodium_calcium_exchanger_current (dimensionless).
    * pcnst[27] is gamma in component sodium_calcium_exchanger_current (dimensionless).
    * pcnst[28] is Km_Ca in component sodium_calcium_exchanger_current (millimolar).
    * pcnst[29] is Km_Nai in component sodium_calcium_exchanger_current (millimolar).
    * pcnst[30] is g_pCa in component calcium_pump_current (picoA_per_picoF).
    * pcnst[31] is K_pCa in component calcium_pump_current (millimolar).
    * pcnst[32] is g_pK in component potassium_pump_current (nanoS_per_picoF).
    * pcnst[33] is k1_prime in component calcium_dynamics (per_millimolar2_per_millisecond).
    * pcnst[34] is k2_prime in component calcium_dynamics (per_millimolar_per_millisecond).
    * pcnst[35] is k3 in component calcium_dynamics (per_millisecond).
    * pcnst[36] is k4 in component calcium_dynamics (per_millisecond).
    * pcnst[37] is EC in component calcium_dynamics (millimolar).
    * pcnst[38] is max_sr in component calcium_dynamics (dimensionless).
    * pcnst[39] is min_sr in component calcium_dynamics (dimensionless).
    * pcnst[40] is V_rel in component calcium_dynamics (per_millisecond).
    * pcnst[41] is V_xfer in component calcium_dynamics (per_millisecond).
    * pcnst[42] is K_up in component calcium_dynamics (millimolar).
    * pcnst[43] is V_leak in component calcium_dynamics (per_millisecond).
    * pcnst[44] is Vmax_up in component calcium_dynamics (millimolar_per_millisecond).
    * pcnst[45] is Buf_c in component calcium_dynamics (millimolar).
    * pcnst[46] is K_buf_c in component calcium_dynamics (millimolar).
    * pcnst[47] is Buf_sr in component calcium_dynamics (millimolar).
    * pcnst[48] is K_buf_sr in component calcium_dynamics (millimolar).
    * pcnst[49] is Buf_ss in component calcium_dynamics (millimolar).
    * pcnst[50] is K_buf_ss in component calcium_dynamics (millimolar).
    * pcnst[51] is V_sr in component calcium_dynamics (micrometre3).
    * pcnst[52] is V_ss in component calcium_dynamics (micrometre3).
    */
   double pcnst[53]; 
   pcnst[0] = 8314.472;
   pcnst[1] = 310;
   pcnst[2] = 96485.3415;
   pcnst[3] = 0.185;
   pcnst[4] = 0.016404;
   pcnst[5] = 10;
   pcnst[6] = 1000;
   pcnst[7] = 1;
   pcnst[8] = 52;
   pcnst[9] = 0.03;
   pcnst[10] = 5.4;
   pcnst[11] = 140;
   pcnst[12] = 2;
   pcnst[13] = 5.405;
   pcnst[14] = 0.153;
   /*
   pcnst[15] = 0.392;  //endo
   pcnst[15] = 0.098;  //mid
   pcnst[15] = 0.392;  //Epi
   */
   pcnst[16] = 14.838;
   pcnst[17] = 0.00029;
   pcnst[18] = 0.0000398;
   pcnst[19] = 0.000592;
   /*
   pcnst[20] = 0.073;   //endo
   pcnst[20] = 0.294;   //mid
   pcnst[20] = 0.294;   //Epi
   */
   pcnst[21] = 2.724;
   pcnst[22] = 1;
   pcnst[23] = 40;
   pcnst[24] = 1000;
   pcnst[25] = 0.1;
   pcnst[26] = 2.5;
   pcnst[27] = 0.35;
   pcnst[28] = 1.38;
   pcnst[29] = 87.5;
   pcnst[30] = 0.1238;
   pcnst[31] = 0.0005;
   pcnst[32] = 0.0146;
   pcnst[33] = 0.15;
   pcnst[34] = 0.045;
   pcnst[35] = 0.06;
   pcnst[36] = 0.005;
   pcnst[37] = 1.5;
   pcnst[38] = 2.5;
   pcnst[39] = 1;
   pcnst[40] = 0.102;
   pcnst[41] = 0.0038;
   pcnst[42] = 0.00025;
   pcnst[43] = 0.00036;
   pcnst[44] = 0.006375;
   pcnst[45] = 0.2;
   pcnst[46] = 0.001;
   pcnst[47] = 10;
   pcnst[48] = 0.3;
   pcnst[49] = 0.4;
   pcnst[50] = 0.00025;
   pcnst[51] = 0.001094;
   pcnst[52] = 0.00005468;
   
   
   c1 = pcnst[2]/(pcnst[0]*pcnst[1]); 
   c2 = pcnst[9]; 
   c3 = -1/c1; 
   c4 = -c3*log(pcnst[11]);
   c5 = -c3*log(pcnst[10]);
   c6 = -c3*log(pcnst[10]+pcnst[9]*pcnst[11]);
   //
   c8 = -0.5*c3*log(pcnst[12]);
   c9 = -pcnst[3]/(pcnst[4]*pcnst[2]);
   
   double g_Ks_Endo_Epi= 0.392; 
   double g_Ks_Mid= 0.098; 
   double g_to_Endo = 0.073;
   double g_to_Mid_Epi = 0.294;
   
   //cA[0] =  c9*g_Ks_Endo_Epi; 
   //cA[1] =  c9*g_Ks_Mid; 
   //cA[2] =  c9*g_Ks_Endo_Epi; 

   g_Ks[0] =  g_Ks_Endo_Epi; 
   g_Ks[1] =  g_Ks_Mid; 
   g_Ks[2] =  g_Ks_Endo_Epi; 
   
   //cB[0] =  c9*g_to_Endo; 
   //cB[1] =  c9*g_to_Mid_Epi; 
   //cB[2] =  c9*g_to_Mid_Epi; 

   g_to[0] =  g_to_Endo; 
   g_to[1] =  g_to_Mid_Epi; 
   g_to[2] =  g_to_Mid_Epi; 
   
   c10= 1/(0.5*c9);
   c7 =  (0.50*pcnst[19]);
   c11= pcnst[14]*sqrt(pcnst[10]/5.4);
   c12= pcnst[13]*sqrt(pcnst[10]/5.4);
   c13= pcnst[3]/(2.0*pcnst[52]*pcnst[2]*c9);
   c14 = pcnst[51]*pcnst[40]/(pcnst[52]*c9); 
   c15 = -pcnst[52]/pcnst[4]; 
   c16 = pcnst[51]/pcnst[4];
   c17 = pcnst[35]/(pcnst[33]*pcnst[34])*c9;
   c18 = pcnst[34]*pcnst[38]/c9 ;
   c19  = -pcnst[34]*(pcnst[38]-pcnst[39])/c9; 
   c20  = pcnst[16]; 
   c21  = pcnst[17]; 
   c22  = 1/c9; 
   c23  = pcnst[41]/(c15*c9); 
   //c24  = pcnst[30]/c10; 
   c24  =  0.5*pcnst[30]; 
   c25  =  1.0/pcnst[23]; 
   c26  =  1.0/pcnst[31]; 
   c27  =  1.0/pcnst[42]; 
   c28  =  1.0/sqrt(pcnst[45]*pcnst[46]); 
   c29  =  pcnst[46]*c28; 
   c30  =  1.0/pcnst[37]; 
   c31  =  1.0/sqrt(pcnst[47]*pcnst[48]); 
   c32  =  pcnst[48]*c31; 
   c33  =  1.0/sqrt(pcnst[49]*pcnst[50]); 
   c34  =  pcnst[50]*c33; 
   c36  =  pcnst[36]/c9; 
   c40  =  pcnst[40]/c9; 
   c43  =  pcnst[43]/c9; 
   c44  =  pcnst[44]/c9; 
   
   
   P_NaK= pcnst[21]; 
   g_NaL =0.0; 
   
   f1 = c1; 
   f2 =  -2.0*pcnst[10]/(pcnst[10]+pcnst[22]);
   f3 =  ((CUBE(pcnst[29])+CUBE(pcnst[11]))*(pcnst[28]+pcnst[12]))/(pcnst[24]*pcnst[12]); 
   f4 =  f3*pcnst[25]; 
   f5 =  pcnst[27]*f1; 
   f6 =  (CUBE(pcnst[11])*pcnst[26]/pcnst[12]);
   f7 = pcnst[18]*pcnst[2]*f1;
   f7a = 0.5*pcnst[18]*pcnst[2];
   f9  = 4.0*pcnst[12];
   f9a = 4.0*pcnst[12]*f7a;
   f10 = pcnst[32];

   cnst.c2 =c2;  
   cnst.c3 =c3;  
   cnst.c4 =c4;  
   cnst.c5 =c5;  
   cnst.c6 =c6;  
   cnst.c7 =c7;  
   cnst.c8 =c8;  
   cnst.c9 =c9;  
   cnst.c11 =c11;  
   cnst.c13 =c13;  
   cnst.c14 =c14;  
   cnst.c15 =c15;  
   cnst.c16 =c16;  
   cnst.c17 =c17;  
   cnst.c18 =c18;  
   cnst.c19 =c19;  
   cnst.c20 =c20;  
   cnst.c21 =c21;  
   cnst.c22 =c22;  
   cnst.c23 =c23;  
   cnst.c24 =c24;  
   cnst.c25 =c25;  
   cnst.c26 =c26;  
   cnst.c27 =c27;  
   cnst.c28 =c28;  
   cnst.c29 =c29;  
   cnst.c30 =c30;  
   cnst.c31 =c31;  
   cnst.c32 =c32;  
   cnst.c33 =c33;  
   cnst.c34 =c34;  
   cnst.c36 =c36;  
   cnst.c40 =c40;  
   cnst.c43 =c43;  
   cnst.c44 =c44;  
   cnst.twenty = 20.0; 
   cnst.eighty = 80.0; 
   cnst.point6 = 0.6;
   cnst.point4 = 0.4;
   set_SP(cnst); 

}

double get_c9() { return cnst.c9; }
double fv0(double Vm, void *p) 
{ 
   double expV1 = exp(-f1*Vm); 
   return f2/(1.0+0.1245*exp(-0.1*Vm*f1)+0.0353*expV1);
}
double fv1(double Vm, void *p) 
{ 
   double expV1 = exp(-f1*Vm); 
   double expV5 = exp(f5*Vm); 
   double expV = expV5*expV1; 
   double fva = 1.0/(f3+ f4*expV);
   return expV5*fva ; 
}
double fv2(double Vm, void *p) 
{ 
   double expV1 = exp(-f1*Vm); 
   double expV5 = exp(f5*Vm); 
   double expV = expV5*expV1; 
   double fva = 1.0/(f3+ f4*expV);
   return expV* fva*f6 ; 
}
double fv3(double Vm, void *p) 
{ 
   double x = 2.0*f1*(Vm-15.0);
   double x2 =x*x;  
   if (x2 < 1e-4) return (1 + 0.5*x + x2/12.0-x2*x2/720)*f7a ;
   return f7a*x/(1.0-exp(-x));
}
double fv4(double Vm, void *p) 
{
   double x = 2*f1*(Vm-15.0);
   double x2 =x*x;  
   if (x2 < 1e-4) return (1 - 0.5*x + x2/12.0-x2*x2/720)*f9a ;
   return f9a*x/(exp(x)-1.0);
}
double fv5(double Vm, void *p) 
{
   return f10/(1.0+exp((25.0 - Vm)/5.98));
}
double fv6(double dV0, void *p) 
{
   double xx  =  (3.0*exp(0.0002*dV0 + 0.02)+exp(0.1*dV0 - 1.0))/(1.0+exp( -0.5*dV0))*(10.0+10*exp(0.06*dV0 -12.0));
   double fv6 = c12/(1.0+xx);
   return fv6; 
   
}

//double ss = cnst.c2*(states[Na_i]/states[K_i]);
//double dV2 = Vm-cnst.c3*log(states[K_i]) -cnst.c3*log(1+ss) -c6;
//double dV2 = Vm-cnst.c3*log(states[K_i])-cnst.c3*logSeries(cnst.c2*states[Na_i]/states[K_i]) -cnst.c6;
//double dV2 = Vm-cnst.c3*log(states[K_i]+cnst.c2*states[Na_i]) -cnst.c6;
//#define logSeries(x) ((x)*(1.0+(x)*(-0.5+(x)/3.0)))

#define logSeries(x)    (log(1+(x)) )

#define fv0_m 10
#define fv1_m  6
#define fv2_m  8
#define fv3_m 11
#define fv4_m  9
#define fv5_m  8
#define fv6_m  7
#define fv5_l  5
#define fv6_l 11

void fv05Special(double v, double *fv)
{
  double fv0_a[]={-1.45499489004044e+00,-2.24272406249374e-03, 2.87019818110585e-05,-2.93525992825876e-07, 1.86878131970215e-09,-6.15332407381053e-12,-3.07155144145422e-14, 2.10687633250732e-15,-1.57414794358693e-17,-1.69089871040937e-19};
  double fv1_a[]={ 1.57571330412697e-04, 2.41325169024597e-06, 1.46053164590980e-08, 6.04272220479134e-11, 2.60192856479526e-13, 6.02170007043617e-16};
  double fv2_a[]={ 5.40504807620401e+02,-1.19531891610964e+01, 1.18653494080701e-01,-6.00926106818024e-04, 4.00772142804151e-07, 2.94210996176821e-08,-2.55280464729594e-10,-1.52370496566743e-12};
  double fv3_a[]={ 1.03965405517718e+00, 4.60580107291175e-02, 7.91075918239993e-04, 4.30241562974461e-06,-4.97022959453108e-08,-6.96639048514165e-10, 2.49915783972090e-12, 7.81666161976579e-14, 1.36854440783922e-16,-3.65847291145083e-18,-1.74506487282187e-20};
  double fv4_a[]={ 2.55684105370893e+01,-7.81445921919197e-01, 6.31701538469039e-03, 3.41107297247441e-05,-3.79918615957248e-07,-5.22361846923760e-09, 1.18578052624143e-11, 4.55066993718199e-13, 1.89136344887477e-15};
  double fv5_a[]={ 2.19831706117241e-04, 2.16983094649991e-05, 9.83170991434406e-07, 2.61347101393462e-08, 4.31098192118690e-10, 4.32295749601499e-12, 2.40194725811093e-14, 5.64722239204646e-17, 
                  1.00000000000000e+00,-6.61573242254081e-02, 2.00764738834124e-03,-2.65106715048650e-05, 2.24306976181154e-07};
  double sum0 = 0; for (int j=fv0_m-1;j>=0    ;j--)sum0  =  fv0_a[j]+v*sum0;
  double sum1 = 0; for (int j=fv1_m-1;j>=0    ;j--)sum1  =  fv1_a[j]+v*sum1;
  double sum2 = 0; for (int j=fv2_m-1;j>=0    ;j--)sum2  =  fv2_a[j]+v*sum2;
  double sum3 = 0; for (int j=fv3_m-1;j>=0    ;j--)sum3  =  fv3_a[j]+v*sum3;
  double sum4 = 0; for (int j=fv4_m-1;j>=0    ;j--)sum4  =  fv4_a[j]+v*sum4;
  double sum5 = 0; for (int j=fv5_m-1;j>=0    ;j--)sum5  =  fv5_a[j]+v*sum5;
  double sum5d =0; for (int j=fv5_l+fv5_m-1;j>=fv5_m;j--)sum5d =  fv5_a[j] + v*sum5d; 

  fv[0] = sum0; 
  fv[1] = sum1; 
  fv[2] = sum2; 
  fv[3] = sum3; 
  fv[4] = sum4; 
  fv[5] = sum5/sum5d; 
      
}
double fv6Special(double dv)
{
  double fv6_a[]={ 2.97942000558175e-01,-8.89837887994794e-03, 7.99339727741565e-03,-1.75996322715896e-04, 1.50747331544023e-06,-5.91201683806020e-09, 8.93564690058337e-12, 
                  1.00000000000000e+00, 2.17226317175310e-01, 2.17628774116084e-02, 1.59235097855721e-03,-7.48970504345503e-05, 3.38737191083109e-06,-5.36894894634127e-08, 8.84744201829959e-10,-4.46027161749137e-12, 2.07853354660018e-14, 1.95808469577737e-16};
  double sum6 = 0; for (int j=fv6_m-1;j>=0    ;j--)sum6  =  fv6_a[j]+dv*sum6;
  double sum6d =0; for (int j=fv6_l+fv6_m-1;j>=fv6_m;j--)sum6d =  fv6_a[j] + dv*sum6d; 
  return sum6/sum6d; 
}

void fv05Pade(double v, double *fv)
{
  double fv0_a[]={-1.45499489004044e+00,-2.24272406249374e-03, 2.87019818110585e-05,-2.93525992825876e-07, 1.86878131970215e-09,-6.15332407381053e-12,-3.07155144145422e-14, 2.10687633250732e-15,-1.57414794358693e-17,-1.69089871040937e-19};
  double fv1_a[]={ 1.57571330412697e-04, 2.41325169024597e-06, 1.46053164590980e-08, 6.04272220479134e-11, 2.60192856479526e-13, 6.02170007043617e-16};
  double fv2_a[]={ 5.40504807620401e+02,-1.19531891610964e+01, 1.18653494080701e-01,-6.00926106818024e-04, 4.00772142804151e-07, 2.94210996176821e-08,-2.55280464729594e-10,-1.52370496566743e-12};
  double fv3_a[]={ 1.03965405517718e+00, 4.60580107291175e-02, 7.91075918239993e-04, 4.30241562974461e-06,-4.97022959453108e-08,-6.96639048514165e-10, 2.49915783972090e-12, 7.81666161976579e-14, 1.36854440783922e-16,-3.65847291145083e-18,-1.74506487282187e-20};
  double fv4_a[]={ 2.55684105370893e+01,-7.81445921919197e-01, 6.31701538469039e-03, 3.41107297247441e-05,-3.79918615957248e-07,-5.22361846923760e-09, 1.18578052624143e-11, 4.55066993718199e-13, 1.89136344887477e-15};
  double fv5_a[]={ 2.19831706117241e-04, 2.16983094649991e-05, 9.83170991434406e-07, 2.61347101393462e-08, 4.31098192118690e-10, 4.32295749601499e-12, 2.40194725811093e-14, 5.64722239204646e-17, 
                  1.00000000000000e+00,-6.61573242254081e-02, 2.00764738834124e-03,-2.65106715048650e-05, 2.24306976181154e-07};
  double sum0 = 0; for (int j=fv0_m-1;j>=0    ;j--)sum0  =  fv0_a[j]+v*sum0;
  double sum1 = 0; for (int j=fv1_m-1;j>=0    ;j--)sum1  =  fv1_a[j]+v*sum1;
  double sum2 = 0; for (int j=fv2_m-1;j>=0    ;j--)sum2  =  fv2_a[j]+v*sum2;
  double sum3 = 0; for (int j=fv3_m-1;j>=0    ;j--)sum3  =  fv3_a[j]+v*sum3;
  double sum4 = 0; for (int j=fv4_m-1;j>=0    ;j--)sum4  =  fv4_a[j]+v*sum4;
  double sum5 = 0; for (int j=fv5_m-1;j>=0    ;j--)sum5  =  fv5_a[j]+v*sum5;
  double sum5d =0; for (int j=fv5_l+fv5_m-1;j>=fv5_m;j--)sum5d =  fv5_a[j] + v*sum5d; 

  fv[0] = sum0; 
  fv[1] = sum1; 
  fv[2] = sum2; 
  fv[3] = sum3; 
  fv[4] = sum4; 
  fv[5] = sum5/sum5d; 
      
}
double fv6Pade(double dv)
{
  double fv6_a[]={ 2.97942000558175e-01,-8.89837887994794e-03, 7.99339727741565e-03,-1.75996322715896e-04, 1.50747331544023e-06,-5.91201683806020e-09, 8.93564690058337e-12, 
                  1.00000000000000e+00, 2.17226317175310e-01, 2.17628774116084e-02, 1.59235097855721e-03,-7.48970504345503e-05, 3.38737191083109e-06,-5.36894894634127e-08, 8.84744201829959e-10,-4.46027161749137e-12, 2.07853354660018e-14, 1.95808469577737e-16};
  double sum6 = 0; for (int j=fv6_m-1;j>=0    ;j--)sum6  =  fv6_a[j]+dv*sum6;
  double sum6d =0; for (int j=fv6_l+fv6_m-1;j>=fv6_m;j--)sum6d =  fv6_a[j] + dv*sum6d; 
  return sum6/sum6d; 
}
/*
void fv05Exact(double Vm, double *fv)
{
   fv[0] = fv0(Vm,NULL); 
   fv[1] = fv1(Vm,NULL); 
   fv[2] = fv2(Vm,NULL); 
   fv[3] = fv3(Vm,NULL); 
   fv[4] = fv4(Vm,NULL); 
   fv[5] = fv5(Vm,NULL); 
}
double fv6Exact(double dv)
{
   return fv6(dv,NULL); 
} */

void update_nonGate(void *ptr, double dt, struct CellTypeParms *cellTypeParms, int nCells, int *cellTypeVector, double *VM, int offset, double **state, double *dVdt)
{
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

  double *__Na_i = state[Na_i]+offset; 
  double *__Ca_i = state[Ca_i]+offset; 
  double *__Ca_ss = state[Ca_ss]+offset; 
  double *__Ca_SR = state[Ca_SR]+offset; 
  double *__fCass = state[fCass]+offset; 
  double *__dVK_i = state[dVK_i]+offset; 
  double *__R_prime = state[R_prime]+offset; 

  //void fv05General(double Vm, double *fv);
  //double fv6General(double dv);
  void (*fv05Func)(double Vm, double *fv);
  double (*fv6Func)(double dv);

  fv05Func = fv05Pade; 
  fv6Func  = fv6Pade; 
  double fv[6];
 int cellType=-1; 
 for (int ii=0;ii<nCells;ii++) 
 {
   double P_NaK,g_Ks,g_to,g_NaL,midK_i,midNa_i,alphaK_i,alphaNa_i,cK_i,cNa_i; 
   if (cellType != cellTypeVector[ii])
   {
   cellType = cellTypeVector[ii]; 
   P_NaK = cellTypeParms[cellType].P_NaK; 
   g_Ks  = cellTypeParms[cellType].g_Ks; 
   g_to  = cellTypeParms[cellType].g_to; 
   g_NaL = cellTypeParms[cellType].g_NaL; 
   midK_i = cellTypeParms[cellType].midK_i; 
   midNa_i = cellTypeParms[cellType].midNa_i; 
   alphaK_i = 1.0/midK_i; 
   alphaNa_i = 1.0/midNa_i; 
   alphaNa_i = alphaNa_i=1.0; 

   cK_i=  +cnst.c3*log(alphaK_i) -cnst.c5;
   cNa_i= +cnst.c3*log(alphaNa_i)-cnst.c4;
   }

   double Vm = VM[ii]; 
  fv05Func(Vm,fv); 
   double _Na_i    = __Na_i[ii]; 
   double _Ca_i    = __Ca_i[ii]; 
   double _Ca_ss   = __Ca_ss[ii]; 
   double _Ca_SR   = __Ca_SR[ii]; 
   double _fCass   = __fCass[ii]; 
   double _dVK_i   = __dVK_i[ii]; 
   double _R_prime = __R_prime[ii]; 

   
   //double *states = cell[ii].state; 

   double dVR = 0.0; 
   double itmp0,itmp5,itmp6 ; 
   {
     //double fv1=fit[1]->afunc(Vm, fit[1]->aparms); 
     //double fv2=fit[2]->afunc(Vm, fit[2]->aparms); 
     double fv1 = fv[1]; 
     double fv2 = fv[2]; 
     double x1 = _Ca_i*cnst.c26; 
     double x2 = SQ(_Ca_i*cnst.c27); 
     double x3 = SQ(_Ca_i*cnst.c28+cnst.c29); 
     double sigm1 = sigm(x1); 
     double sigm2 = sigm(x2); 
     double sigm3 = sigm(x3); 
     double dV3 = Vm- 0.5*cnst.c3*log(_Ca_i) - cnst.c8;

     itmp0 = (CUBE(_Na_i)*fv1-_Ca_i*fv2); 
     double itmp4 = (cnst.c7*dV3+cnst.c24*sigm1); 
     itmp5=  (cnst.c43*(_Ca_i -  _Ca_SR)+cnst.c44*sigm2);      
     itmp6 = (cnst.c23*(_Ca_ss - _Ca_i));
     __Ca_i[ii]   = _Ca_i + (dt*cnst.c9)*(sigm3*(itmp4-itmp0+itmp6*cnst.c15-itmp5*cnst.c16));
     dVR  -= 2*itmp4;
    }

   {
     double _K_i  = cnst.c9*(_dVK_i-Vm);
     double x0 = _Na_i*cnst.c25; 
     double sigm0 = P_NaK*sigm(x0); 

     double dV0 = Vm -cnst.c3*log(_K_i) -cnst.c5;
     double dV1 = Vm -cnst.c3*log(_Na_i)-cnst.c4;
     double dV2 = dV0-cnst.c3*logSeries(cnst.c2*_Na_i/_K_i)+cnst.c5-cnst.c6; // Assumption:  cnst.c2*_Na_i/_K_i is small; 

     //double fv0=fit[0]->afunc(Vm, fit[0]->aparms); 
     //double fv5=fit[5]->afunc(Vm, fit[5]->aparms); 
     //double fv6=fit[6]->afunc(dV0,fit[6]->aparms); 
     double fv0 = fv[0]; 
     double fv5 = fv[5]; 
     double fv6 = fv6Func(dV0);
     double fd =  fv5  +  fv6; 

     double tmp0 =  (fd +  g_to*rGate[ii]*sGate[ii]+ cnst.c11*Xr1Gate[ii]*Xr2Gate[ii] );
     double tmp1 =  (cnst.c20*CUBE(mGate[ii])*hGate[ii]*jGate[ii]+cnst.c21);
     double tmp2 =  g_Ks*SQ(XsGate[ii]);
     double itmpA = sigm0 * fv0;                          //Sigm0
     double itmp2 = itmp0 - 1.5*itmpA+tmp1*dV1; 
     double itmp3 = itmpA + tmp0*dV0 +tmp2*dV2; 
     double iNaL = g_NaL*CUBE(mGate[ii])*jLGate[ii]*dV1;

     _dVK_i += dt*itmp3;
     dVR    +=  iNaL - itmp2-itmp3;
     __Na_i[ii]  =   _Na_i  +  (dt*cnst.c9)*(iNaL*cnst.c22+itmp2+2.0*itmp0);
   }
   {
     //double fv3=fit[3]->afunc(Vm, fit[3]->aparms); 
     //double fv4=fit[4]->afunc(Vm, fit[4]->aparms); 
     double fv3 = fv[3]; 
     double fv4 = fv[4]; 
     double x4 = SQ(_Ca_SR*cnst.c30); 
     double x5 = SQ(_Ca_SR*cnst.c31+cnst.c32); 
     double x6 = SQ(_Ca_ss*cnst.c33+cnst.c34); 
     double sigm4 = sigm(x4); 
     double sigm5 = sigm(x5); 
     double sigm6 = sigm(x6); 
     double tmp8  = (cnst.c18+cnst.c19*sigm4); //Sigm4
     double tmp9  = tmp8*_Ca_ss+cnst.c36; 
     double itmp1 = dGate[ii]*fGate[ii]*f2Gate[ii]*_fCass*(fv4-_Ca_ss*fv3);
     double sigmCa_ss =   SQ(_Ca_ss)/(tmp8*cnst.c17 + SQ(_Ca_ss));
     double itmp7 = sigmCa_ss*_R_prime*(_Ca_SR - _Ca_ss);

     double t1 = 1.0/(1.0+SQ(20*_Ca_ss)); 
     double mhu = 0.600000*t1+0.4000000;
     double tauR =    1.0/(80.0*t1+2.0);

     __Ca_ss[ii]   = _Ca_ss   + (dt*cnst.c9)*sigm6*(itmp6+itmp7*cnst.c14+itmp1*cnst.c13);    
     __Ca_SR[ii]   = _Ca_SR   + (dt*cnst.c9)*sigm5*(itmp5-cnst.c40*itmp7);
     __R_prime[ii] = _R_prime + (dt*cnst.c9)*(cnst.c36 - tmp9*_R_prime);
     __fCass[ii]   = _fCass   + dt*(mhu - _fCass)*tauR; 
     dVR += itmp1; 
   }
   __dVK_i[ii] = _dVK_i + dt*dVR ; 
   dVdt[ii]  = dVR;
   }
}
