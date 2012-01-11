#include "TT06Func.hh"
#include <stdio.h>
#include <math.h>
#include <cassert>
#include "gsl.h"
#include "TT06Tau.hh"
#include "pade.hh"
#include "mpiUtils.h"

#define C(i) (gsl_vector_get(c,(i)))
#define sigm(x)   ((x)/(1.0+(x)) )

#if (0) 
typedef double (*GATE)(double Vm, double *gate, int cellType, double *A, double *B); 
typedef struct {double (*update)(double Vm, double *x0, int cellType, double *A, double *B); int cellType;  }GateParms;  
double update_fCassGate(double Ca_aa, double *gate, int cellType, double *a, double *b);
double update_Xr1Gate(double Vm, double *gate, int cellType, double *a, double *b);
double update_Xr2Gate(double Vm, double *gate,  int cellType, double *a, double *b);
double update_XsGate(double Vm, double *gate,  int cellType, double *a, double *b);
double update_mGate(double Vm, double *gate, int cellType, double *a, double *b);
double update_hGate(double Vm, double *gate, int cellType, double *a, double *b);
double update_hGateMod(double Vm, double *gate, int cellType, double *a, double *b);
double update_jGate(double Vm, double *gate, int cellType, double *a, double *b);
double update_jGateMod(double Vm, double *gate, int cellType, double *a, double *b);
double update_rGate(double Vm, double *gate, int cellType, double *a, double *b);
double update_dGate(double Vm, double *gate, int cellType, double *a, double *b);
double update_fGate(double Vm, double *gate, int cellType, double *a, double *b);
double update_f2Gate(double Vm, double *gate, int cellType, double *a, double *b);
double update_sGate(double Vm, double *gate, int cellType,double *a, double *b);
static    GATE gate[] ={ update_fCassGate, update_Xr1Gate, update_Xr2Gate,update_XsGate,update_mGate,update_hGate,update_jGate,update_rGate,update_dGate,update_fGate,update_f2Gate, update_sGate}; 
static    char *nameA[] ={ "fCassGateA", "Xr1GateA", "Xr2GateA","XsGateA","mGateA","hGateA","jGateA","rGateA","dGateA","fGateA","f2GateA","sGateA"}; 
static    char *nameB[] ={ "fCassGateB", "Xr1GateB", "Xr2GateB","XsGateB","mGateB","hGateB","jGateB","rGateB","dGateB","fGateB","f2GateB","sGateA"}; 
#endif 

double Xr1Mhu(double Vm, void *parms) ;
double Xr2Mhu(double Vm, void *parms) ;
double XsMhu(double Vm, void *parms) ;
double mMhu(double Vm, void *parms) ;
double hjMhu(double Vm, void *parms) ;
double rMhu(double Vm, void *parms) ;
double dMhu(double Vm, void *parms) ;
double fMhu(double Vm, void *parms) ;
double f2Mhu(double Vm, void *parms) ;
double sMhu(double Vm, void *parms) ;
double fCassMhu(double Vm, void *parms);

double Xr1TauR(double Vm, void *parms) ;
double Xr2TauR(double Vm, void *parms) ;
double XsTauR(double Vm, void *parms) ;
double mTauR(double Vm, void *parms) ;
double hTauR(double Vm, void *parms) ;
double jTauR(double Vm, void *parms) ;
double rTauR(double Vm, void *parms) ;
double dTauR(double Vm, void *parms) ;
double fTauR(double Vm, void *parms) ;
double f2TauR(double Vm, void *parms) ;
double sTauR0(double Vm, void *parms) ;
double sTauR1(double Vm, void *parms) ;
double fCassTauR(double Vm, void *parms) ;

static double tauRMax[] = { 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
static int cnt = 0; 
static TauRecipParms *jParms;
static TauRecipParms *hParms;
static PADE *fit[35]; 
static double cA[3],cB[3]; 
static double c1,c2,c3,c4,c5,c6,c7,c8,c9;
static double c10,c11,c12,c13,c14,c15,c16,c17,c18,c19;
static double c20,c21,c22,c23,c24,c25,c26,c27,c28,c29;
static double c30,c31,c32,c33,c34,c36,c40,c43,c44;
static double f1,f2,f3,f4,f5,f6,f7,f7a,f9,f9a,f10; 
typedef double (*OVF)(double x, void *parms) ;
static    OVF  mhuFunc[] ={ fCassMhu, mMhu, hjMhu, hjMhu ,Xr1Mhu, Xr2Mhu,XsMhu,rMhu,dMhu,fMhu,f2Mhu, sMhu}; 
static    char *mhuName[] ={ "fCassMhu", "mMhu", "hjMhu","hjMhu","Xr1Mhu", "Xr2Mhu","XsMhu","rMhu","dMhu","fMhu","f2Mhu", "sMhu"}; 
static    OVF  tauRFunc[]  ={  fCassTauR ,  mTauR , hTauR ,  jTauR,  Xr1TauR,   Xr2TauR,  XsTauR,  rTauR,  dTauR,  fTauR,  f2TauR,   sTauR0,   sTauR1}; 
static    char *tauRName[] ={ "fCassTauR", "mTauR", "hTauR", "jTauR","Xr1TauR", "Xr2TauR","XsTauR","rTauR","dTauR","fTauR","f2TauR", "sTauR0", "sTauR1"}; 
void initState(double *STATES,int cellType)
{
if (cellType == 0) 
{
//STATES[Vmembrane] = -86.709;
STATES[K_i] = 138.4;
STATES[Na_i] = 10.355;
STATES[Ca_i] = 0.00013;
STATES[Ca_ss] = 0.00036;
STATES[Ca_SR] = 3.715;
STATES[R_prime] = 0.9068;
STATES[fCass_gate] = 0.9953;
STATES[Xr1_gate] = 0.00448;
STATES[Xr2_gate] = 0.476;
STATES[Xs_gate] = 0.0087;
STATES[m_gate] = 0.00155;
STATES[h_gate] = 0.7573;
STATES[j_gate] = 0.7225;
STATES[r_gate] = 2.235e-8;
STATES[d_gate] = 3.164e-5;
STATES[f_gate] = 0.8009;
STATES[f2_gate] = 0.9778;
STATES[s_gate] = 0.3212;
}

if (cellType == 1) 
{
//STATES[Vmembrane] = -85.423;
STATES[K_i] = 138.52;
STATES[Na_i] = 10.132;
STATES[Ca_i] = 0.000153;
STATES[Ca_ss] = 0.00042;
STATES[Ca_SR] = 4.272;
STATES[R_prime] = 0.8978;
STATES[fCass_gate] = 0.9942;
STATES[Xr1_gate] = 0.0165;
STATES[Xr2_gate] = 0.473;
STATES[Xs_gate] = 0.0174;
STATES[m_gate] = 0.00165;
STATES[h_gate] = 0.749;
STATES[j_gate] = 0.6788;
STATES[r_gate] = 2.347e-8;
STATES[d_gate] = 3.288e-5;
STATES[f_gate] = 0.7026;
STATES[f2_gate] = 0.9526;
STATES[s_gate] = 0.999998;
}

if (cellType==2) 
{
//STATES[Vmembrane] = -85.23;
STATES[K_i] = 136.89;
STATES[Na_i] = 8.604;
STATES[Ca_i] = 0.000126;
STATES[Ca_ss] = 0.00036;
STATES[Ca_SR] = 3.64;
STATES[R_prime] = 0.9073;
STATES[fCass_gate] = 0.9953;
STATES[Xr1_gate] = 0.00621;
STATES[Xr2_gate] = 0.4712;
STATES[Xs_gate] = 0.0095;
STATES[m_gate] = 0.00172;
STATES[h_gate] = 0.7444;
STATES[j_gate] = 0.7045;
STATES[r_gate] = 2.42e-8;
STATES[d_gate] = 3.373e-5;
STATES[f_gate] = 0.7888;
STATES[f2_gate] = 0.9755;
STATES[s_gate] = 0.999998;
}
}
void initCnst()
{
   /*
    * cnst[0] is R in component membrane (joule_per_mole_kelvin).
    * cnst[1] is T in component membrane (kelvin).
    * cnst[2] is F in component membrane (coulomb_per_millimole).
    * cnst[3] is Cm in component membrane (microF).
    * cnst[4] is V_c in component membrane (micrometre3).
    * cnst[5] is stim_start in component membrane (millisecond).
    * cnst[6] is stim_period in component membrane (millisecond).
    * cnst[7] is stim_duration in component membrane (millisecond).
    * cnst[8] is stim_amplitude in component membrane (picoA_per_picoF).
    * cnst[9] is P_kna in component reversal_potentials (dimensionless).
    * cnst[10] is K_o in component potassium_dynamics (millimolar).
    * cnst[11] is Na_o in component sodium_dynamics (millimolar).
    * cnst[12] is Ca_o in component calcium_dynamics (millimolar).
    * cnst[13] is g_K1 in component inward_rectifier_potassium_current (nanoS_per_picoF).
    * cnst[14] is g_Kr in component rapid_time_dependent_potassium_current (nanoS_per_picoF).
    * cnst[15] is g_Ks in component slow_time_dependent_potassium_current (nanoS_per_picoF).
    * cnst[16] is g_Na in component fast_sodium_current (nanoS_per_picoF).
    * cnst[17] is g_bna in component sodium_background_current (nanoS_per_picoF).
    * cnst[18] is g_CaL in component L_type_Ca_current (litre_per_farad_second).
    * cnst[19] is g_bca in component calcium_background_current (nanoS_per_picoF).
    * cnst[20] is g_to in component transient_outward_current (nanoS_per_picoF).
    * cnst[21] is P_NaK in component sodium_potassium_pump_current (picoA_per_picoF).
    * cnst[22] is K_mk in component sodium_potassium_pump_current (millimolar).
    * cnst[23] is K_mNa in component sodium_potassium_pump_current (millimolar).
    * cnst[24] is K_NaCa in component sodium_calcium_exchanger_current (picoA_per_picoF).
    * cnst[25] is K_sat in component sodium_calcium_exchanger_current (dimensionless).
    * cnst[26] is alpha in component sodium_calcium_exchanger_current (dimensionless).
    * cnst[27] is gamma in component sodium_calcium_exchanger_current (dimensionless).
    * cnst[28] is Km_Ca in component sodium_calcium_exchanger_current (millimolar).
    * cnst[29] is Km_Nai in component sodium_calcium_exchanger_current (millimolar).
    * cnst[30] is g_pCa in component calcium_pump_current (picoA_per_picoF).
    * cnst[31] is K_pCa in component calcium_pump_current (millimolar).
    * cnst[32] is g_pK in component potassium_pump_current (nanoS_per_picoF).
    * cnst[33] is k1_prime in component calcium_dynamics (per_millimolar2_per_millisecond).
    * cnst[34] is k2_prime in component calcium_dynamics (per_millimolar_per_millisecond).
    * cnst[35] is k3 in component calcium_dynamics (per_millisecond).
    * cnst[36] is k4 in component calcium_dynamics (per_millisecond).
    * cnst[37] is EC in component calcium_dynamics (millimolar).
    * cnst[38] is max_sr in component calcium_dynamics (dimensionless).
    * cnst[39] is min_sr in component calcium_dynamics (dimensionless).
    * cnst[40] is V_rel in component calcium_dynamics (per_millisecond).
    * cnst[41] is V_xfer in component calcium_dynamics (per_millisecond).
    * cnst[42] is K_up in component calcium_dynamics (millimolar).
    * cnst[43] is V_leak in component calcium_dynamics (per_millisecond).
    * cnst[44] is Vmax_up in component calcium_dynamics (millimolar_per_millisecond).
    * cnst[45] is Buf_c in component calcium_dynamics (millimolar).
    * cnst[46] is K_buf_c in component calcium_dynamics (millimolar).
    * cnst[47] is Buf_sr in component calcium_dynamics (millimolar).
    * cnst[48] is K_buf_sr in component calcium_dynamics (millimolar).
    * cnst[49] is Buf_ss in component calcium_dynamics (millimolar).
    * cnst[50] is K_buf_ss in component calcium_dynamics (millimolar).
    * cnst[51] is V_sr in component calcium_dynamics (micrometre3).
    * cnst[52] is V_ss in component calcium_dynamics (micrometre3).
    */
   double cnst[53]; 
   cnst[0] = 8314.472;
   cnst[1] = 310;
   cnst[2] = 96485.3415;
   cnst[3] = 0.185;
   cnst[4] = 0.016404;
   cnst[5] = 10;
   cnst[6] = 1000;
   cnst[7] = 1;
   cnst[8] = 52;
   cnst[9] = 0.03;
   cnst[10] = 5.4;
   cnst[11] = 140;
   cnst[12] = 2;
   cnst[13] = 5.405;
   cnst[14] = 0.153;
   /*
   cnst[15] = 0.392;  //endo
   cnst[15] = 0.098;  //mid
   cnst[15] = 0.392;  //Epi
   */
   cnst[16] = 14.838;
   cnst[17] = 0.00029;
   cnst[18] = 0.0000398;
   cnst[19] = 0.000592;
   /*
   cnst[20] = 0.073;   //endo
   cnst[20] = 0.294;   //mid
   cnst[20] = 0.294;   //Epi
   */
   cnst[21] = 2.724;
   cnst[22] = 1;
   cnst[23] = 40;
   cnst[24] = 1000;
   cnst[25] = 0.1;
   cnst[26] = 2.5;
   cnst[27] = 0.35;
   cnst[28] = 1.38;
   cnst[29] = 87.5;
   cnst[30] = 0.1238;
   cnst[31] = 0.0005;
   cnst[32] = 0.0146;
   cnst[33] = 0.15;
   cnst[34] = 0.045;
   cnst[35] = 0.06;
   cnst[36] = 0.005;
   cnst[37] = 1.5;
   cnst[38] = 2.5;
   cnst[39] = 1;
   cnst[40] = 0.102;
   cnst[41] = 0.0038;
   cnst[42] = 0.00025;
   cnst[43] = 0.00036;
   cnst[44] = 0.006375;
   cnst[45] = 0.2;
   cnst[46] = 0.001;
   cnst[47] = 10;
   cnst[48] = 0.3;
   cnst[49] = 0.4;
   cnst[50] = 0.00025;
   cnst[51] = 0.001094;
   cnst[52] = 0.00005468;
   
   
   c1 = cnst[2]/(cnst[0]*cnst[1]); 
   c2 = cnst[9]; 
   c3 = -1/c1; 
   c4 = -c3*log(cnst[11]);
   c5 = -c3*log(cnst[10]);
   c6 = -c3*log(cnst[10]+cnst[9]*cnst[11]);
   //
   c8 = -0.5*c3*log(cnst[12]);
   c9 = -cnst[3]/(cnst[4]*cnst[2]);
   
   double g_Ks_Endo_Epi= 0.392; 
   double g_Ks_Mid= 0.098; 
   double g_to_Endo = 0.073;
   double g_to_Mid_Epi = 0.294;
   
   cA[0] =  c9*g_Ks_Endo_Epi; 
   cA[1] =  c9*g_Ks_Mid; 
   cA[2] =  c9*g_Ks_Endo_Epi; 
   
   cB[0] =  c9*g_to_Endo; 
   cB[1] =  c9*g_to_Mid_Epi; 
   cB[2] =  c9*g_to_Mid_Epi; 
   
   c10= 1/(0.5*c9);
   c7 =  (0.50*cnst[19]*c9);
   c11= c9*cnst[14]*sqrt(cnst[10]/5.4);
   c12= cnst[13]*sqrt(cnst[10]/5.4);
   c13= cnst[3]/(2.0*cnst[52]*cnst[2]);
   c14 = cnst[51]*cnst[40]/cnst[52]; 
   c15 = -cnst[52]/cnst[4]; 
   c16 = cnst[51]/cnst[4];
   c17 = cnst[35]/(cnst[33]*cnst[34]);
   c18 = cnst[34]*cnst[38] ;
   c19  = -cnst[34]*(cnst[38]-cnst[39]); 
   c20  = c9*cnst[16]; 
   c21  = c9*cnst[17]; 
   c22  = -1/c9; 
   c23  = cnst[41]/c15; 
   c24  = cnst[30]/c10; 
   c25  =  1.0/cnst[23]; 
   c26  =  1.0/cnst[31]; 
   c27  =  1.0/cnst[42]; 
   c28  =  1.0/sqrt(cnst[45]*cnst[46]); 
   c29  =  cnst[46]*c28; 
   c30  =  1.0/cnst[37]; 
   c31  =  1.0/sqrt(cnst[47]*cnst[48]); 
   c32  =  cnst[48]*c31; 
   c33  =  1.0/sqrt(cnst[49]*cnst[50]); 
   c34  =  cnst[50]*c33; 
   c36  =  cnst[36]; 
   c40  =  cnst[40]; 
   c43  =  cnst[43]; 
   c44  =  cnst[44]; 
   
   f1 = c1; 
   f2 =  -2.0*c9*cnst[21]*cnst[10]/(cnst[10]+cnst[22]);
   f3 =  ((CUBE(cnst[29])+CUBE(cnst[11]))*(cnst[28]+cnst[12]))/(cnst[24]*cnst[12]*c9); 
   f4 =  f3*cnst[25]; 
   f5 =  cnst[27]*f1; 
   f6 =  (CUBE(cnst[11])*cnst[26]/cnst[12]);
   f7 = cnst[18]*cnst[2]*f1;
   f7a = 0.5*cnst[18]*cnst[2];
   f9 = 4.0*cnst[12];
   f9a = 4.0*cnst[12]*f7a;
   f10 = c9*cnst[32];
   jParms =makeTauRecipParms(-48.85,-17.6,jTauRecip); 
   hParms =makeTauRecipParms(-64.20,-23.3,hTauRecip); 
   
}
double get_c9() { return c9; }
double fv0(double Vm, void *parm) 
{ 
   double expV1 = exp(-f1*Vm); 
   return f2/(1.0+0.1245*exp(-0.1*Vm*f1)+0.0353*expV1);
}
double fv1(double Vm, void *parm) 
{ 
   double expV1 = exp(-f1*Vm); 
   double expV5 = exp(f5*Vm); 
   double expV = expV5*expV1; 
   double fva = 1.0/(f3+ f4*expV);
   return expV5*fva ; 
}
double fv2(double Vm, void *parm) 
{ 
   double expV1 = exp(-f1*Vm); 
   double expV5 = exp(f5*Vm); 
   double expV = expV5*expV1; 
   double fva = 1.0/(f3+ f4*expV);
   return expV* fva*f6 ; 
}
double fv3(double Vm, void *parm) 
{ 
   double x = 2.0*f1*(Vm-15.0);
   double x2 =x*x;  
   if (x2 < 1e-4) return (1 + 0.5*x + x2/12.0-x2*x2/720)*f7a ;
   return f7a*x/(1.0-exp(-x));
}
double fv4(double Vm, void *parm) 
{
   double x = 2*f1*(Vm-15.0);
   double x2 =x*x;  
   if (x2 < 1e-4) return (1 - 0.5*x + x2/12.0-x2*x2/720)*f9a ;
   return f9a*x/(exp(x)-1.0);
}
double fv5(double Vm, void *parm) 
{
   return f10/(1.0+exp((25.0 - Vm)/5.98));
}
double fv6(double dV0, void *parm) 
{
   double xx  =  (3.0*exp(0.0002*dV0 + 0.02)+exp(0.1*dV0 - 1.0))/(1.0+exp( -0.5*dV0))*(10.0+10*exp(0.06*dV0 -12.0));
   double fv6 = c9*c12/(1.0+xx);

   //static FILE *file=NULL; 
   //static double t =0.0; 
   //if (file==NULL) {file =fopen("fv6.data","w");  t=0.0; }
    //fprintf(file,"%f %e %e %e\n",t,dV0,xx,fv6); 
    //t+=0.02; 
   return fv6; 
   
}
void makeFit(double tol, double V0, double V1, double deltaV , int maxOrder, int maxCost,int mod) 
{
   int nPade=0;; 
   fit[nPade++]=padeApprox("fv0",fv0,NULL,0,tol,deltaV,V0,V1,maxOrder,maxCost); 
   fit[nPade++]=padeApprox("fv1",fv1,NULL,0,tol,deltaV,V0,V1,maxOrder,maxCost); 
   fit[nPade++]=padeApprox("fv2",fv2,NULL,0,tol,deltaV,V0,V1,maxOrder,maxCost); 
   fit[nPade++]=padeApprox("fv3",fv3,NULL,0,tol,deltaV,V0,V1,maxOrder,maxCost); 
   fit[nPade++]=padeApprox("fv4",fv4,NULL,0,tol,deltaV,V0,V1,maxOrder,maxCost); 
   fit[nPade++]=padeApprox("fv5",fv5,NULL,0,tol,deltaV,V0,V1,maxOrder,maxCost); 
   fit[nPade++]=padeApprox("fv6",fv6,NULL,0,tol,deltaV,0,130,maxOrder,maxCost); 

   int k=0; 
   fit[nPade++]=padeApprox( mhuName[k], mhuFunc[k],NULL,0,tol,deltaV,0.0,2.0,maxOrder,maxCost); 
   fit[nPade++]=padeApprox(tauRName[k],tauRFunc[k],NULL,0,tol,deltaV,0.0,2.0,maxOrder,maxCost); 
   for (k=1;k<12;k++) 
   {
     fit[nPade++]=padeApprox( mhuName[k], mhuFunc[k],NULL,0,tol,deltaV,V0,V1,maxOrder,maxCost); 
     fit[nPade++]=padeApprox(tauRName[k],tauRFunc[k],NULL,0,tol,deltaV,V0,V1,maxOrder,maxCost); 
   }
   fit[nPade++]=padeApprox(tauRName[k],tauRFunc[k],NULL,0,tol,deltaV,V0,V1,maxOrder,maxCost); 
/*
    int k=0; 
   parms.update = gate[k] ;
   fit[nPade++]=padeApprox(mhuName[k],getGateA,(void*)&parms,sizeof(parms),tol,deltaV,0.0,2.0,maxOrder,maxCost); 
   fit[nPade++]=padeApprox(tauRName[k],getGateB,(void*)&parms,sizeof(parms),tol,deltaV,0.0,2.0,maxOrder,maxCost); 
   for (k=1;k<12;k++) 
   {
      parms.update = gate[k] ;
      fit[nPade++]=padeApprox(mhuName[k],getGateA,(void*)&parms,sizeof(parms),tol,deltaV,V0,V1,maxOrder,maxCost); 
      fit[nPade++]=padeApprox(tauRName[k],getGateB,(void*)&parms,sizeof(parms),tol,deltaV,V0,V1,maxOrder,maxCost); 
   }
   parms.cellType = 1; 
   parms.update = gate[11]; 
   fit[nPade++]=padeApprox("sGateA_1",getGateA,(void*)&parms,sizeof(parms),tol,deltaV,V0,V1,maxOrder,maxCost); 
   fit[nPade++]=padeApprox("sGateB_1",getGateB,(void*)&parms,sizeof(parms),tol,deltaV,V0,V1,maxOrder,maxCost); 

*/
   if (getRank(0)==0) 
   {
      FILE *file=fopen("pade.data","w"); 
      for (int  i=0;i<nPade;i++) 
      {
         padeErrorInfo(*fit[i]); 
         padeWrite(file,*fit[i]); 
      }
      fclose(file); 
   }
}

double defaultVoltage(int cellType)
{
   switch (cellType)
   {
     case 0:
      return -86.709;
      break;
     case 1:
      return -85.423;
      break;
     case 2:
      return -85.23;
      break;
     default:
      assert(false);
   }
}



#define logSeries(x) ((x)*(1.0+(x)*(-0.5+(x)/3.0)))
double computeUpdates(double dt, double Vm, double* STATES, int cellType, double *rates)
{

   
   double x[8]; 
   x[0] = STATES[Na_i]*c25; 
   x[1] = STATES[Ca_i]*c26; 
   x[2] = SQ(STATES[Ca_i]*c27); 
   x[3] = SQ(STATES[Ca_i]*c28+c29); 
   x[4] = SQ(STATES[Ca_SR]*c30); 
   x[5] = SQ(STATES[Ca_SR]*c31+c32); 
   x[6] = SQ(STATES[Ca_ss]*c33+c34); 
   double sigm0 = sigm(x[0]); 
   double sigm1 = sigm(x[1]); 
   double sigm2 = sigm(x[2]); 
   double sigm3 = sigm(x[3]); 
   double sigm4 = sigm(x[4]); 
   double sigm5 = sigm(x[5]); 
   double sigm6 = sigm(x[6]); 
   
   double dV0 = Vm -c3*log(STATES[K_i]) -c5;
   double dV1 = Vm -c3*log(STATES[Na_i])-c4;
   double dV2 = dV0-c3*logSeries(c2*STATES[Na_i]/STATES[K_i])+c5-c6; // Assumption:  c2*STATES[Na_i]/STATES[K_i] is small; 
   double dV3 = Vm- 0.5*c3*log(STATES[Ca_i]) - c8;
   //double ss = c2*(STATES[Na_i]/STATES[K_i]);
   //double dV2 =Vm-c3*log(STATES[K_i]) -c3*log(1+ss) -c6;
   //double dV2 =Vm-c3*log(STATES[K_i])-c3*logSeries(c2*STATES[Na_i]/STATES[K_i]) -c6;
   //double dV2 = 1.0*Vm -c3*log(STATES[K_i]+c2*STATES[Na_i]) -c6;
   
   
   
   double fv[5]; 
   for (int i=0;i<5;i++) fv[i]=fit[i]->afunc(Vm, fit[i]->aparms); 
   double fd = fit[5]->afunc(Vm, fit[5]->aparms) +fit[6]->afunc(dV0, fit[6]->aparms); 

   double tmp0 =    fd +  cB[cellType]*STATES[r_gate]*STATES[s_gate]+ c11*STATES[Xr1_gate]*STATES[Xr2_gate] ;
   double tmp1 =  c20*CUBE(STATES[m_gate])*STATES[h_gate]*STATES[j_gate]+c21;
   double tmp2 =  cA[cellType]*SQ(STATES[Xs_gate]);
   double tmp8  = c18+c19*sigm4; //Sigm4
   
   double sigmCa_ss =   SQ(STATES[Ca_ss])/(tmp8*c17 + SQ(STATES[Ca_ss]));
   double  tmp9 =    tmp8*STATES[Ca_ss]+c36; 
   
   double itmpA =  sigm0 * fv[0];                          //Sigm0
   double itmp0 = CUBE(STATES[Na_i])*fv[1]-STATES[Ca_i]*fv[2]; 
   double itmp1 = STATES[d_gate]*STATES[f_gate]*STATES[f2_gate]*STATES[fCass_gate]*(fv[4]-STATES[Ca_ss]*fv[3]);
   double itmp2 = itmp0 - 1.5*itmpA+tmp1*dV1; 
   double itmp3 = itmpA + tmp0*dV0 +tmp2*dV2; 
   double itmp4 = c7*dV3+c24*sigm1; 

   double itmp5=  c43*(STATES[Ca_i] - STATES[Ca_SR])+c44*sigm2;      
   double itmp6 = c23*(STATES[Ca_ss] - STATES[Ca_i]);
   double itmp7 = sigmCa_ss*STATES[R_prime]*(STATES[Ca_SR] - STATES[Ca_ss]);
   
   double dVdt  = itmp1+itmp2*c22+itmp3*c22-itmp4*c10;

   double Ca_SS     = STATES[Ca_ss]; 
   rates[K_i]     = (itmp3);
   rates[Na_i]    = (itmp2+2.0*itmp0);
   rates[Ca_i]    = (sigm3*(itmp4-itmp0+itmp6*c15-itmp5*c16));
   rates[Ca_ss]   = (sigm6*(itmp6+itmp7*c14+itmp1*c13));    
   rates[Ca_SR]   = (sigm5*(itmp5-c40*itmp7));
   rates[R_prime] = (c36 - tmp9*STATES[R_prime]);


// Update Gates; 
   
   double mhu,tauR; 
   double *gate = STATES+gateOffset; 
   double *grate = rates+gateOffset; 
   PADE **gatefit = fit + 7; 
   int i=0;
   mhu=gatefit[2*i]  ->afunc(Ca_SS,gatefit[2*i  ]->aparms);    //fCassGate
   tauR=gatefit[2*i+1]->afunc(Ca_SS,gatefit[2*i+1]->aparms); 
   grate[i] =  (mhu - gate[i])*tauR; 
   if ( tauR > tauRMax[i]) tauRMax[i] = tauR; 
   i++; 

   mhu=gatefit[2*i]  ->afunc(Vm,gatefit[2*i  ]->aparms);       //mGate
   tauR=gatefit[2*i+1]->afunc(Vm,gatefit[2*i+1]->aparms); 
   grate[i] =  (mhu - gate[i])*tauR; 
   grate[i]  *= (1-exp(-dt*tauR))/(dt*tauR);                   //mtauR can be very large   ~1140.0 use Rush_Larson to integrate. 
   if ( tauR > tauRMax[i]) tauRMax[i] = tauR; 
   i++;

   mhu=gatefit[2*i]  ->afunc(Vm,gatefit[2*i  ]->aparms);      // hGate 
   tauR=gatefit[2*i+1]->afunc(Vm,gatefit[2*i+1]->aparms); 
   grate[i] =  (mhu - gate[i])*tauR; 
   if ( tauR > tauRMax[i]) tauRMax[i] = tauR; 
   i++;

   //mhu=gatefit[2*i]  ->afunc(Vm,gatefit[2*i  ]->aparms);    //jGate        note hMhu and  hMhu are equal. 
   tauR=gatefit[2*i+1]->afunc(Vm,gatefit[2*i+1]->aparms); 
   grate[i] =  (mhu - gate[i])*tauR; 
   if ( tauR > tauRMax[i]) tauRMax[i] = tauR; 

   for (i=4;i<11;i++)                    // other Gates. 
   {
      mhu=gatefit[2*i]  ->afunc(Vm,gatefit[2*i  ]->aparms); 
      tauR=gatefit[2*i+1]->afunc(Vm,gatefit[2*i+1]->aparms); 
      grate[i] =  (mhu - gate[i])*tauR; 
      if ( tauR > tauRMax[i]) tauRMax[i] = tauR; 
   }
   mhu=gatefit[2*i]  ->afunc(Vm,gatefit[2*i  ]->aparms);                //sGate
   tauR=gatefit[2*i+1+(cellType != 0)  ]->afunc(Vm,gatefit[2*i+1+(cellType != 0)  ]->aparms);  //Note  sTauR depends on celltype
   grate[i] =  (mhu - gate[i])*tauR; 
   if ( tauR > tauRMax[i]) tauRMax[i] = tauR; 

   return dVdt; 
}
void write_tauRMax()
{
	for (int i=0;i<12;i++) fprintf(stderr,"%d %s %14.8f\n",i,tauRName[i],tauRMax[i]); 
}
double fCassMhu(double Ca_SS, void *parms) 
{ 
   double t1 = 1.0/(1.0+SQ(20*Ca_SS)); 
   double mhu = 0.600000*t1+0.4000000;
   return mhu;
}
double fCassTauR(double Ca_SS, void *parms) 
{ 
   double t1 = 1.0/(1.0+SQ(20*Ca_SS)); 
   double tauR =    1.0/(80.0*t1+2.0);
   return tauR;
}
double Xr1Mhu(double Vm, void *parms) 
{ 
   double mhu=1.0/(1.0+(exp(((-26.0 - Vm)/7.0))));
   return mhu ; 
}
double Xr1TauR(double Vm, void *parms) 
{ 
   double t1 = (1.0+(exp(((-45.0 - Vm)/10.0))))/450;
   double t2 = (1.0+(exp(((Vm+30.0)/11.5000))))/6.0;
   double tauR =  t1*t2;
   return tauR;
}
double Xr2Mhu(double Vm, void *parms) 
{ 
   double mhu=1.0/(1.0+(exp(((Vm+88.0)/24.0))));
   return mhu;
}
double Xr2TauR(double Vm, void *parms) 
{ 
   double t1 = (1.0+(exp(((-60.0 - Vm)/20.0))))/3.0;
   double t2 = (1.0+(exp(((Vm - 60.0)/20.0))))/1.120;
   double tauR =  t1*t2;
   return tauR;
}
double XsMhu(double Vm, void *parms) 
{ 
   double mhu = 1.0/(1.0+(exp(((-5.0 - Vm)/14.0))));
   return mhu;
}
double XsTauR(double Vm, void *parms) 
{ 
   double t1 = 1400.00/ sqrt(1.0+exp((5.0 - Vm)/6.0));
   double t2 = 1.0/(1.0+(exp(((Vm - 35.0)/15.0))));
   double tauR  =  1.0/(t1*t2+80.0);
   return tauR;
}
double mMhu(double Vm, void *parms) 
{ 
   double mhu = 1.0/SQ(1.0+exp((-56.8600 - Vm)/9.03000));
   return mhu;
}
double mTauR(double Vm, void *parms) 
{ 
   double t1  = 1.0/(1.0+(exp(((- 60.0 - Vm)/5.0))));
   double t2  =  0.10000/(1.0+(exp(((Vm+35.0)/5.0))))+0.100000/(1.0+(exp(((Vm - 50.0)/200.0))));
   double tauR =  1.0/(t1*t2);
   return tauR;
}
double hjMhu(double Vm, void *parms) 
{ 
   double mhu = 1.0/SQ((1.0+(exp(((Vm+71.5500)/7.43000)))));
   return mhu;
}
double hTauR(double Vm, void *parms) 
{ 
   double t1  = (Vm<- 40.0 ?  0.0570000*(exp((- (Vm+80.0)/6.80000))) : 0.0);
   double t2  = (Vm<- 40.0 ?  2.70000*(exp(( 0.0790000*Vm)))+ 310000.*(exp(( 0.348500*Vm))) : 0.770000/( 0.130000*(1.0+(exp(((Vm+10.6600)/- 11.1000))))));
   double tauR = (t1+t2);
   return tauR;
}
double jTauR(double Vm, void *parms) 
{ 
   double t1  = (Vm < -40.0 ? (( ( - 25428.0*(exp(( 0.244400*Vm))) -  6.94800e-06*(exp(( - 0.0439100*Vm))))*(Vm+37.7800))/1.0)/(1.0+(exp(( 0.311000*(Vm+79.2300))))) : 0.0);
   double t2 = (Vm < -40.0 ? ( 0.0242400*(exp(( - 0.0105200*Vm))))/(1.0+(exp(( - 0.137800*(Vm+40.1400))))) : ( 0.600000*(exp(( 0.0570000*Vm))))/(1.0+(exp(( - 0.100000*(Vm+32.0))))));
   double tauR  = (t1+t2);
   return tauR;
}
double rMhu(double Vm, void *parms) 
{ 
   double mhu = 1.0/(1.0+(exp(((20.0 - Vm)/6.0))));
   return mhu;
}
double rTauR(double Vm, void *parms) 
{ 
   double tau =  9.50000*(exp((- SQ((Vm+40.0)))/1800.00))+0.800000;
   double tauR = 1.0/tau; 
   return tauR;
}
double dMhu(double Vm, void *parms) 
{ 
   double mhu = 1.0/(1.0+(exp(((- 8.0 - Vm)/7.50000))));
   return mhu;
}
double dTauR(double Vm, void *parms) 
{ 
   double t1  = 1.40000/(1.0+(exp(((- 35.0 - Vm)/13.0))))+0.250000;
   double t2 = 1.40000/(1.0+(exp(((Vm+5.0)/5.0))));
   double t3 = 1.0/(1.0+(exp(((50.0 - Vm)/20.0))));
   double tauR =  1/(t1*t2+t3);
   return tauR;
}
double fMhu(double Vm, void *parms) 
{ 
   double mhu = 1.0/(1.0+(exp(((Vm+20.0)/7.0))));
   return mhu;
}
double fTauR(double Vm, void *parms) 
{ 
   double tau =  1102.50*(exp((- SQ(Vm+27.0)/225.0)))+200.0/(1.0+(exp(((13.0 - Vm)/10.0))))+180.0/(1.0+(exp(((Vm+30.0)/10.0))))+20.0;
   double tauR = 1/tau; 
   return tauR;
}
double f2Mhu(double Vm, void *parms) 
{ 
   double mhu = 0.670000/(1.0+(exp(((Vm+35.0)/7.0))))+0.330000;
   return mhu;
}
double f2TauR(double Vm, void *parms) 
{ 
   double tau =  562.0*exp(-SQ((Vm+27.0))/240.0)+31.0/(1.0+(exp(((25.0 - Vm)/10.0))))+80.0/(1.0+(exp(((Vm+30.0)/10.0))));
   double tauR = 1/tau; 
   return tauR;
}
double sMhu(double Vm, void *parms) 
{ 
   double mhu = 1.0/(1.0+(exp(((Vm+28.0)/5.0))));
   return mhu;
}
double sTauR0(double Vm, void *parms) 
{ 
   double   tau =  1000.0*(exp((-SQ(Vm+67.0)/1000.0)))+8.0;
   double tauR = 1/tau; 
   return tauR;
}
double sTauR1(double Vm, void *parms) 
{ 
   double   tau =  85.0*(exp((- SQ(Vm+45.0)/320.0)))+5.0/(1.0+(exp(((Vm - 20.0)/5.0))))+3.0;
   double tauR = 1/tau; 
   return tauR;
}
#if (0) 
double getGateA(double Vm, void *parms)
{
   double x0,A,B;
   GateParms *p = (GateParms *)parms; 
   x0 = 0; p->update(Vm, &x0, p->cellType, &A, &B);
   return A; 
}
double getGateB(double Vm, void *parms)
{
   double x0,A,B;
   GateParms *p = (GateParms *)parms; 
   x0 = 0; p->update(Vm, &x0, p->cellType, &A, &B);
   return B; 
}
double update_Xr1Gate(double Vm, double *gate, int cellType, double *a, double *b)
{
   double mhu = 1.0/(1.0+(exp(((- 26.0 - Vm)/7.0))));
   double t1 = 450.0/(1.0+(exp(((- 45.0 - Vm)/10.0))));
   double t2 = 6.0/(1.0+(exp(((Vm+30.0)/11.5000))));
   double tau =  t1*t2;
   double rate = (mhu - *gate)/tau;
   *a =  mhu; 
   *b =  1.0/tau; 
   return rate; 
}
double update_Xr2Gate(double Vm, double *gate, int cellType,  double *a, double *b)
{
   double mhu = 1.0/(1.0+(exp(((Vm+88.0)/24.0))));
   double t1 = 3.0/(1.0+(exp(((- 60.0 - Vm)/20.0))));
   double t2 = 1.12000/(1.0+(exp(((Vm - 60.0)/20.0))));
   double tau =  t1*t2;
   double rate = (mhu - *gate)/tau;
   *a =  mhu; 
   *b =  1.0/tau; 
   return rate; 
}
double update_XsGate(double Vm, double *gate, int cellType,  double *a, double *b)
{
   double mhu = 1.0/(1.0+(exp(((- 5.0 - Vm)/14.0))));
   double t1 = 1400.00/ sqrt(1.0+exp((5.0 - Vm)/6.0));
   double t2 = 1.0/(1.0+(exp(((Vm - 35.0)/15.0))));
   double tauR  =  1.0/(t1*t2+80.0);
   *a =  mhu; 
   *b =  tauR; 
   double rate = (mhu - *gate)*tauR;
   return rate; 
}
double update_mGate(double Vm, double *gate, int cellType, double *a, double *b)
{
   double mhu = 1.0/SQ(1.0+exp((-56.8600 - Vm)/9.03000));
   double t1  = 1.0/(1.0+(exp(((- 60.0 - Vm)/5.0))));
   double t2  =  0.10000/(1.0+(exp(((Vm+35.0)/5.0))))+0.100000/(1.0+(exp(((Vm - 50.0)/200.0))));
   double tauR =  1.0/(t1*t2);
   *a = mhu; 
   *b = tauR ;
   double rate = (mhu-*gate)*tauR;
   return rate; 
}
double update_hGate(double Vm, double *gate, int cellType, double *a, double *b)
{
   double mhu = 1.0/SQ((1.0+(exp(((Vm+71.5500)/7.43000)))));
   double t1  = (Vm<- 40.0 ?  0.0570000*(exp((- (Vm+80.0)/6.80000))) : 0.0);
   double t2  = (Vm<- 40.0 ?  2.70000*(exp(( 0.0790000*Vm)))+ 310000.*(exp(( 0.348500*Vm))) : 0.770000/( 0.130000*(1.0+(exp(((Vm+10.6600)/- 11.1000))))));
   double tauR = (t1+t2);
   *a =  mhu; 
   *b =  tauR; 
   double rate = (mhu - *gate)*tauR;
   return rate; 
}
double update_hGateMod(double Vm, double *gate, int cellType, double *a, double *b)
{
   double mhu = 1.0/SQ((1.0+(exp(((Vm+71.5500)/7.43000)))));
   double dtauR,ddtauR; 
   double tauR = TauRecipMod(Vm,hParms,&dtauR,&ddtauR); 
   *a =  mhu; 
   *b =  tauR; 
   double rate = (mhu - *gate)*tauR;
   return rate; 
}
double update_jGate(double Vm, double *gate, int cellType, double *a, double *b)
{
   double mhu = 1.0/SQ((1.0+(exp(((Vm+71.5500)/7.43000)))));
   double t1  = (Vm < -40.0 ? (( ( - 25428.0*(exp(( 0.244400*Vm))) -  6.94800e-06*(exp(( - 0.0439100*Vm))))*(Vm+37.7800))/1.0)/(1.0+(exp(( 0.311000*(Vm+79.2300))))) : 0.0);
   double t2 = (Vm < -40.0 ? ( 0.0242400*(exp(( - 0.0105200*Vm))))/(1.0+(exp(( - 0.137800*(Vm+40.1400))))) : ( 0.600000*(exp(( 0.0570000*Vm))))/(1.0+(exp(( - 0.100000*(Vm+32.0))))));
   double tauR  = (t1+t2);
   *a =  mhu; 
   *b =  tauR; 
   double rate = (mhu - *gate)*tauR;
   return rate; 
}
double update_jGateMod(double Vm, double *gate, int cellType, double *a, double *b)
{
   double mhu = 1.0/SQ((1.0+(exp(((Vm+71.5500)/7.43000)))));
   double dtauR,ddtauR; 
   double tauR = TauRecipMod(Vm,jParms,&dtauR,&ddtauR); 
   *a =  mhu; 
   *b =  tauR; 
   double rate = (mhu - *gate)*tauR;
   return rate; 
}
double update_rGate(double Vm, double *gate, int cellType, double *a, double *b)
{
   double mhu = 1.0/(1.0+(exp(((20.0 - Vm)/6.0))));
   double tau =  9.50000*(exp((- SQ((Vm+40.0)))/1800.00))+0.800000;
   double tauR = 1.0/tau; 
   *a =  mhu; 
   *b =  tauR; 
   double rate = (mhu - *gate)*tauR;
   return rate; 
}
double update_dGate(double Vm, double *gate, int cellType, double *a, double *b)
{
   double mhu = 1.0/(1.0+(exp(((- 8.0 - Vm)/7.50000))));
   double t1  = 1.40000/(1.0+(exp(((- 35.0 - Vm)/13.0))))+0.250000;
   double t2 = 1.40000/(1.0+(exp(((Vm+5.0)/5.0))));
   double t3 = 1.0/(1.0+(exp(((50.0 - Vm)/20.0))));
   double tauR =  1/(t1*t2+t3);
   *a =  mhu; 
   *b =  tauR; 
   double rate = (mhu - *gate)*tauR;
   return rate; 
}
double update_fGate(double Vm, double *gate, int cellType, double *a, double *b)
{
   double mhu = 1.0/(1.0+(exp(((Vm+20.0)/7.0))));
   double tau =  1102.50*(exp((- SQ(Vm+27.0)/225.0)))+200.0/(1.0+(exp(((13.0 - Vm)/10.0))))+180.0/(1.0+(exp(((Vm+30.0)/10.0))))+20.0;
   double tauR = 1/tau; 
   *a =  mhu; 
   *b =  tauR; 
   double rate = (mhu - *gate)*tauR;
   return rate; 
}
double update_f2Gate(double Vm, double *gate, int cellType, double *a, double *b)
{
   double mhu = 0.670000/(1.0+(exp(((Vm+35.0)/7.0))))+0.330000;
   double tau =  562.0*exp(-SQ((Vm+27.0))/240.0)+31.0/(1.0+(exp(((25.0 - Vm)/10.0))))+80.0/(1.0+(exp(((Vm+30.0)/10.0))));
   double tauR = 1/tau; 
   *a =  mhu; 
   *b =  tauR; 
   double rate = (mhu - *gate)*tauR;
   return rate; 
}
double update_sGate(double Vm, double *gate, int cellType,double *a, double *b)
{

   double rate,mhu,tau;
   if (cellType==0) 
   {
      mhu = 1.0/(1.0+(exp(((Vm+28.0)/5.0))));
      tau =  1000.0*(exp((- (pow((Vm+67.0), 2.0))/1000.0)))+8.0;
   }
   else
   {
      mhu = 1.0/(1.0+(exp(((Vm+20.0)/5.0))));
      tau =  85.0*(exp((- (pow((Vm+45.0), 2.0))/320.0)))+5.0/(1.0+(exp(((Vm - 20.0)/5.0))))+3.0;
   }
   double tauR = 1/tau; 
   *a =  mhu; 
   *b =  tauR; 
   rate = (mhu - *gate)*tauR;
   return rate; 
}
double update_fCassGate(double Ca_SS, double *gate, int cellType, double *a, double *b)
{
   double t1 = 1.0/(1.0+SQ(20*Ca_SS)); 
   double mhu = 0.600000*t1+0.4000000;
   double tauR =    1.0/(80.0*t1+2.0);
   *a =  mhu; 
   *b =  tauR; 
   double rate = (mhu - *gate)*tauR;
   //*gate += rate*dt;
    return rate; 
/*
   static FILE *file=NULL; 
   static double t =0.0; 
   if (file==NULL) {file =fopen("fCassGate.data","w");  t=0.0; }
    fprintf(file,"%f %e %e %e %e %e\n",t,Ca_SS,mhu,tau,*a,*b); 
   t+=dt;
*/
}
#endif 
