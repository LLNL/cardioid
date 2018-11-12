#include <math.h>
#include <stdio.h>
#include <omp.h> 
#include "TT06Func.h" 
#include "TT06NonGates.h" 
#define sigm(x)   ((x)/(1.0+(x)) )


static struct nonGateCnst cnst;
static  double f1,f2,f3,f4,f5,f6,f7,f7a,f9,f9a,f10,c12; 
void initNonGateCnst()
{
 double g_Ks[3],g_Kr[3],g_to[3]; 
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
   

   g_Ks[0] =  g_Ks_Endo_Epi; 
   g_Ks[1] =  g_Ks_Mid; 
   g_Ks[2] =  g_Ks_Endo_Epi; 
   g_Kr[0] = 0.153;
   g_Kr[1] = 0.153;
   g_Kr[2] = 0.153;
   

   g_to[0] =  g_to_Endo; 
   g_to[1] =  g_to_Mid_Epi; 
   g_to[2] =  g_to_Mid_Epi; 
   
   c7 =  pcnst[19];
   //c11= pcsnt[14]*sqrt(pcnst[10]/5.4);
   c11= sqrt(pcnst[10]/5.4);
   c12= pcnst[13]*sqrt(pcnst[10]/5.4);
   c13= pcnst[3]/(2.0*pcnst[52]*pcnst[2]*c9);
   c14 = pcnst[51]/pcnst[52]; 
   c15 = -pcnst[52]/pcnst[4]; 
   c16 = pcnst[51]/pcnst[4];
   c17 = pcnst[35]/(pcnst[33]*pcnst[34])*c9;
   c18 = pcnst[34]*pcnst[38]/c9 ;
   c19  = -pcnst[34]*(pcnst[38]-pcnst[39])/c9; 
   c20  = pcnst[16]; 
   c21  = pcnst[17]; 
   c22  = 1/c9; 
   c23  = pcnst[41]/(c15*c9); 
   c24  =  pcnst[30]; 
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
void  update_fCassGate(double dt, int nCells, double *Ca_ss, double *Ca_i, double *g, double *mhu_a, double *tauR_a)
{

   for (int ii=0;ii<nCells;ii++)
   {
      double t1 = 1.0/(1.0+SQ(20*Ca_ss[ii]));
      double mhu = 0.600000*t1+0.4000000;
      double tauR =    1.0/(80.0*t1+2.0);
      g[ii] +=  (mhu - g[ii])*tauR*dt;
   }
}
void  update_fCassGateRice(double dt, int nCells, double *Ca_ss, double *Ca_i, double *g, double *mhu_a, double *tauR_a)
{
   for (int ii=0;ii<nCells;ii++)
   {
      double xCa_ss = 1000*Ca_ss[ii];
      double xCa_i  = 1000*Ca_i[ii] ;
      double mhu    = 0.6/(1.0 + SQ(xCa_ss)/sqrt(xCa_ss) + SQ(xCa_i)/sqrt(xCa_i)) + 0.4;
      double tauR   = 0.005/mhu;
      g[ii] +=  (mhu - g[ii])*tauR*dt;
   }
}
void sampleLog(struct LogParms *logParms, int nCells, int offset, int *cellTypeVector, double *VM, double **state)
{

  if (logParms->loop % 1000 > 0) return; 
  if (logParms->file == NULL) 
  {   
      char filename[64]; 
      int ompID = omp_get_thread_num(); 
      sprintf(filename,"sampleLog%2.2d.data",ompID); 
      logParms->file = fopen(filename,"w"); 
      logParms->loop = 0; 
  }
  int ii=offset; 
  double Vm = VM[ii]; 
  double _Na_i = state[Na_i][ii]; 
  double _Ca_i = state[Ca_i][ii]; 
  double _dVK_i = state[dVK_i][ii]; 
  double _K_i  = cnst.c9*(_dVK_i-Vm);
  double mCa=_Ca_i; 
  double mNa=_Na_i; 
  double mK =_K_i; 
  double MCa=_Ca_i; 
  double MNa=_Na_i; 
  double MK =_K_i; 
  
  for (ii=offset;ii<nCells+offset;ii++) 
  {
       Vm = VM[ii]; 
       _Na_i = state[Na_i][ii]; 
       _Ca_i = state[Ca_i][ii]; 
       _dVK_i = state[dVK_i][ii]; 
       _K_i  = cnst.c9*(_dVK_i-Vm);
       if (mCa < _Ca_i) mCa=_Ca_i; 
       if (mNa < _Na_i) mNa=_Na_i; 
       if (mK < _K_i) mK=_K_i; 
       if (MCa > _Ca_i) MCa=_Ca_i; 
       if (MNa > _Na_i) MNa=_Na_i; 
       if (MK >  _K_i ) MK =_K_i; 
  }
  fprintf(logParms->file,"%d %e %e %e ",logParms->loop,mCa,mK,mNa); 
  fprintf(logParms->file,"%e %e %e\n",MCa,MK,MNa); 
  fflush(logParms->file); 
  logParms->loop++; 
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

void update_nonGate(void *fit, double dt, struct CellTypeParms *cellTypeParms, int nCells, double *VM, int offset, double **state, double *dVdt)
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

  void fv05General(void *fit,  double Vm, double *fv);
  double fv6General(void *fit, double dv);
  void (*fv05Func)(void *fit, double Vm, double *fv);
  double (*fv6Func)(void *fit, double dv);

  fv05Func = fv05General; 
  fv6Func  = fv6General; 
  double fv[6];
 double  c_K1   =  1; 
 double  c_Na   =  1 *cnst.c20; 
 double  c_bNa  =  1*cnst.c21; 
 double  c_CaL  =  1; 
 double  c_bCa  =  1*cnst.c7; 
 double  c_NaCa =  1; 
 double  c_pCa  =  1*cnst.c24; 
 double  c_pK   =  1; 
 double  c_up   = 1*cnst.c44;
 double  c_leak = 1*cnst.c43;
 double  c_xfer = 1*cnst.c23;
 double  c_rel  = 1  *cnst.c40;
 //CURRENT_SCALES *cS = currentScales; 
 //printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",cS->K1,cS->Na,cS->bNa,cS->CaL,cS->bCa,cS->NaCa,cS->pCa,cS->pK,cS->NaK,cS->Ks,cS->Kr,cS->to,cS->NaL,cS->up,cS->leak,cS->xfer,cS->rel); 
 double c_NaK,c_Ks,c_Kr,c_to,c_NaL;

 double midK_i,midNa_i,alphaK_i,alphaNa_i,cK_i,cNa_i; 
 c_NaK = 1 *cellTypeParms->P_NaK; 
 c_Ks  = 1 *cellTypeParms->g_Ks; 
 c_Kr  = 1 *cellTypeParms->g_Kr; 
 c_to  = 1 *cellTypeParms->g_to; 
 c_NaL = 1 *cellTypeParms->g_NaL; 
 
 midK_i = cellTypeParms->midK_i; 
 midNa_i = cellTypeParms->midNa_i; 
 alphaK_i = 1.0/midK_i; 
 alphaNa_i = 1.0/midNa_i; 
 cK_i=  +cnst.c3*log(alphaK_i) -cnst.c5;
 cNa_i= +cnst.c3*log(alphaNa_i)-cnst.c4;

 for (int ii=0;ii<nCells;ii++) 
 {

   double Vm = VM[ii]; 
   fv05Func(fit, Vm,fv); 
   double _Na_i    = __Na_i[ii]; 
   double _Ca_i    = __Ca_i[ii]; 
   double _Ca_ss   = __Ca_ss[ii]; 
   double _Ca_SR   = __Ca_SR[ii]; 
   double _fCass   = __fCass[ii]; 
   double _dVK_i   = __dVK_i[ii]; 
   double _R_prime = __R_prime[ii]; 

   
   //double *states = cell[ii].state; 

   double dVR = 0.0; 
   //double itmp0,itmp5,itmp6 ; 
   double I_K1,I_Kr,I_Ks,I_Na,I_bNa,I_CaL,I_bCa,I_to,I_NaK,I_NaCa,I_pCa,I_pK,I_NaL,I_leak,I_up,I_rel,I_xfer; 
   double I_sum,I_delta;   
//  Update Ca concentration; 
   {
     double fv1 = fv[1]; 
     double fv2 = fv[2]; 
     double x1 = _Ca_i*cnst.c26; 
     double x2 = SQ(_Ca_i*cnst.c27); 
     double x3 = SQ(_Ca_i*cnst.c28+cnst.c29); 
     double sigm1 = sigm(x1); 
     double sigm2 = sigm(x2); 
     double sigm3 = sigm(x3); 
     double dV3 = Vm- 0.5*cnst.c3*log(_Ca_i) - cnst.c8;

     I_NaCa = c_NaCa*(CUBE(_Na_i)*fv1-_Ca_i*fv2); 

     I_bCa = c_bCa*dV3;
     I_pCa = c_pCa*sigm1;    
     I_up  =  cnst.c44*sigm2;        
     I_leak=  cnst.c43*(_Ca_SR - _Ca_i);
     I_xfer = cnst.c23*(_Ca_ss - _Ca_i);
     I_delta = I_leak-I_up;    // I_detal = -itmp5  
     I_sum =   I_bCa+I_pCa; // itmp4 =I_bCa+I_pCa; 
     __Ca_i[ii]   = _Ca_i + (dt*cnst.c9)*(sigm3*(0.5*I_sum-I_NaCa+I_xfer*cnst.c15+I_delta*cnst.c16));
     dVR  -= I_sum;
     //if (ii %4 ==0) printf("\n%d dVR=%14.12f ", ii,dVR); 
     //else printf("%14.12f ", ii,I_sum); 
    }

   double iK; 
//  Update K and Na concentrations; 
   {
     double _K_i  = cnst.c9*(_dVK_i-Vm);
     double x0 = _Na_i*cnst.c25; 

     double dV0 = Vm -cnst.c3*log(_K_i) -cnst.c5;
     double dV1 = Vm -cnst.c3*log(_Na_i)-cnst.c4;
     double dV2 = dV0-cnst.c3*logSeries(cnst.c2*_Na_i/_K_i)+cnst.c5-cnst.c6; // Assumption:  cnst.c2*_Na_i/_K_i is small; 

     double fv0 = fv[0]; 
     double fv5 = fv[5]; 
     double fv6 = fv6Func(fit,dV0);

     I_NaK = c_NaK*sigm(x0) * fv0;                          // renamed itmpA to I_NaK
     I_pK  = c_pK*fv5*dV0;  
     I_K1  = c_K1*fv6*dV0; 
     I_to  = c_to*rGate[ii]*sGate[ii]*dV0 ;
     I_Kr  = c_Kr*Xr1Gate[ii]*Xr2Gate[ii]*dV0;
     I_Na  = c_Na*CUBE(mGate[ii])*hGate[ii]*jGate[ii]*dV1; 
     I_NaL = c_NaL*CUBE(mGate[ii])*jLGate[ii]*dV1;
     I_bNa = c_bNa*dV1;
     I_Ks  = c_Ks*SQ(XsGate[ii])*dV2;
     //if (ii == 0) printf("c_Ks = %15.12f %15.12f %15.12f %15.12f\n",c_Ks,SQ(XsGate[ii]),dV2, I_Ks); 

     double iNa =  3*I_NaCa - 1.5*I_NaK + I_Na + I_bNa + I_NaL;
     iK =   I_Ks   +     I_NaK + I_pK + I_K1  +  I_to + I_Kr; 


     _dVK_i += dt*iK;
     dVR    -=  iNa + iK - 2.0*I_NaCa;
     __Na_i[ii]  =   _Na_i  +  (dt*cnst.c9)*iNa;
   }
     //if (ii %4 ==0) printf("\n%d dVR=%14.12f ", ii,dVR); else printf("%14.12f ", ii,dVR); 
//if (ii == 0) printf("nonSIMD: dVR= %15.12f I= %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f\n", 
            //dVR, I_NaL, I_NaK, I_NaCa, I_Na, I_bNa, I_pK, I_K1, I_to, I_Kr, I_Ks);

//  Update Ca_SS, Ca_SR, R_prime concentrations and fCass gate; 
   {
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

     I_CaL = c_CaL*dGate[ii]*fGate[ii]*f2Gate[ii]*_fCass*(fv4-_Ca_ss*fv3);  // renamed itmp1 to I_CaL

     double O =   SQ(_Ca_ss)*_R_prime/(tmp8*cnst.c17 + SQ(_Ca_ss));
     I_rel =cnst.c40* O*(_Ca_SR - _Ca_ss);

     __Ca_ss[ii]   = _Ca_ss   + (dt*cnst.c9)*sigm6*(I_xfer+I_rel*cnst.c14+I_CaL*cnst.c13);    
     __Ca_SR[ii]   = _Ca_SR   - (dt*cnst.c9)*sigm5*(I_delta+I_rel);
     __R_prime[ii] = _R_prime + (dt*cnst.c9)*(cnst.c36 - tmp9*_R_prime);

#if fCassForm == TT06
     double t1 = 1.0/(1.0+SQ(20*_Ca_ss)); 
     double mhu = 0.600000*t1+0.4000000;
     double tauR =    1.0/(80.0*t1+2.0);
#endif

#if  fCassForm == RICE 
      double xCa_ss = 1000*_Ca_ss;
      double xCa_i  = 1000*_Ca_i ;
      double mhu    = 0.6/(1.0 + xCa_ss*sqrt(xCa_ss) + xCa_i*sqrt(xCa_i)) + 0.4;
      double tauR   = 0.005/mhu;
#endif
     __fCass[ii]   = _fCass   + dt*(mhu - _fCass)*tauR; 

     dVR += I_CaL; 
   }
//  update voltages 
   __dVK_i[ii] = _dVK_i + dt*dVR ; 
   dVdt[ii]  = dVR;
   //printf("%d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",ii,I_K1,I_Na,I_bNa,I_CaL,I_bCa,I_NaCa,I_pCa,I_pK,I_NaK,I_Ks,I_Kr,I_to,I_NaL,I_leak,I_up,I_rel,I_xfer); 
   }
   //printf("\n"); 
}

