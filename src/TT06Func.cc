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

double Xr1Mhu(double Vm, void *parms) ;
double Xr2Mhu(double Vm, void *parms) ;
double XsMhu(double Vm, void *parms) ;
double mMhu(double Vm, void *parms) ;
double hjMhu(double Vm, void *parms) ;
double rMhu(double Vm, void *parms) ;
double dMhu(double Vm, void *parms) ;
double fMhu(double Vm, void *parms) ;
double f2Mhu(double Vm, void *parms) ;
double sMhu0(double Vm, void *parms) ;
double sMhu1(double Vm, void *parms) ;

double Xr1TauR(double Vm, void *parms) ;
double Xr2TauR(double Vm, void *parms) ;
double XsTauR(double Vm, void *parms) ;
double mTauR(double Vm, void *parms) ;
double hTauR(double Vm, void *parms) ;
double hTauRMod(double Vm, void *parms) ;
double jTauR(double Vm, void *parms) ;
double jTauRMod(double Vm, void *parms) ;
double rTauR(double Vm, void *parms) ;
double dTauR(double Vm, void *parms) ;
double fTauR(double Vm, void *parms) ;
double f2TauR(double Vm, void *parms) ;
double sTauR0(double Vm, void *parms) ;
double sTauR1(double Vm, void *parms) ;

static TauRecipParms *jParms;
static TauRecipParms *hParms;
static PADE *fit[36]; 
static double cA[3],cB[3]; 
static double c1,c2,c3,c4,c5,c6,c7,c8,c9;
static double c10,c11,c12,c13,c14,c15,c16,c17,c18,c19;
static double c20,c21,c22,c23,c24,c25,c26,c27,c28,c29;
static double c30,c31,c32,c33,c34,c36,c40,c43,c44;
static double f1,f2,f3,f4,f5,f6,f7,f7a,f9,f9a,f10; 
typedef double (*OVF)(double x, void *parms) ;
static    OVF  mhuFunc[] ={    mMhu,   hjMhu,  hjMhu , Xr1Mhu, Xr2Mhu,XsMhu,rMhu,dMhu,fMhu,f2Mhu, sMhu0, sMhu1}; 
static  const  char *mhuName[] ={ "mMhu", "hjMhu","hjMhu","Xr1Mhu", "Xr2Mhu","XsMhu","rMhu","dMhu","fMhu","f2Mhu", "sMhu0","sMhu1"}; 
static    OVF  tauRFunc[]  ={ mTauR , hTauR ,  jTauR,  Xr1TauR,   Xr2TauR,  XsTauR,  rTauR,  dTauR,  fTauR,  f2TauR,   sTauR0,   sTauR1}; 
static  const   char *tauRName[] ={ "mTauR", "hTauR", "jTauR","Xr1TauR", "Xr2TauR","XsTauR","rTauR","dTauR","fTauR","f2TauR", "sTauR0", "sTauR1"}; 
void initState(double *STATES,int cellType)
{
if (cellType == 0) 
{
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
   return fv6; 
   
}
void makeFit(double tol, double V0, double V1, double deltaV , int maxOrder, int maxCost,int mod) 
{
   if (mod)
   {
      
        tauRFunc[1]      = hTauRMod;
        tauRName[1]  = "hTauRMod";
        tauRFunc[2]      = jTauRMod;
        tauRName[2]  = "jTauRMod";
   }

   int nPade=0;; 
   int nT = maxOrder+1; 
   if (1) 
   {
   fit[nPade++]=padeApprox("fv0",fv0,NULL,0,tol,deltaV,V0,V1,nT,nT,maxCost); 
   fit[nPade++]=padeApprox("fv1",fv1,NULL,0,tol,deltaV,V0,V1,nT,nT,maxCost); 
   fit[nPade++]=padeApprox("fv2",fv2,NULL,0,tol,deltaV,V0,V1,nT,nT,maxCost); 
   fit[nPade++]=padeApprox("fv3",fv3,NULL,0,tol,deltaV,V0,V1,nT,nT,maxCost); 
   fit[nPade++]=padeApprox("fv4",fv4,NULL,0,tol,deltaV,V0,V1,nT,nT,maxCost); 
   fit[nPade++]=padeApprox("fv5",fv5,NULL,0,tol,deltaV,V0,V1,nT,nT,maxCost); 
   }
   else
   {
   fit[nPade++]=padeApprox("fv0",fv0,NULL,0,tol,deltaV,V0,V1,-1,-6,maxCost); 
   fit[nPade++]=padeApprox("fv1",fv1,NULL,0,tol,deltaV,V0,V1,-1,-6,maxCost); 
   fit[nPade++]=padeApprox("fv2",fv2,NULL,0,tol,deltaV,V0,V1,-1,-6,maxCost); 
   fit[nPade++]=padeApprox("fv3",fv3,NULL,0,tol,deltaV,V0,V1,-1,-6,maxCost); 
   fit[nPade++]=padeApprox("fv4",fv4,NULL,0,tol,deltaV,V0,V1,-1,-6,maxCost); 
   fit[nPade++]=padeApprox("fv5",fv5,NULL,0,tol,deltaV,V0,V1,-3,-6,maxCost); 
   }

   fit[nPade++]=padeApprox("fv6",fv6,NULL,0,tol,deltaV,0,130,nT,nT,maxCost); 

   int k;
   for (k=0;k<11;k++) 
   {
     fit[nPade++]=padeApprox( mhuName[k], mhuFunc[k],NULL,0,tol,deltaV,V0,V1,nT,nT,maxCost); 
     fit[nPade++]=padeApprox(tauRName[k],tauRFunc[k],NULL,0,tol,deltaV,V0,V1,nT,nT,maxCost); 
   }
   fit[nPade++]=padeApprox( mhuName[k], mhuFunc[k],NULL,0,tol,deltaV,V0,V1,nT,nT,maxCost); 
   fit[nPade++]=padeApprox(tauRName[k],tauRFunc[k],NULL,0,tol,deltaV,V0,V1,nT,nT,maxCost); 
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
   double voltage = 0;
   switch (cellType)
   {
     case 0:
      voltage = -86.709;
      break;
     case 1:
      voltage = -85.423;
      break;
     case 2:
      voltage = -85.23;
      break;
     default:
      assert(false);
   }
   return voltage;
}



     //double ss = c2*(states[Na_i]/states[K_i]);
     //double dV2 =Vm-c3*log(states[K_i]) -c3*log(1+ss) -c6;
     //double dV2 =Vm-c3*log(states[K_i])-c3*logSeries(c2*states[Na_i]/states[K_i]) -c6;
     //double dV2 = Vm -c3*log(states[K_i]+c2*states[Na_i]) -c6;

#define logSeries(x) ((x)*(1.0+(x)*(-0.5+(x)/3.0)))
void computeNonGateRates(double dt, int nCells, const double *VM, TT06DevState* States, double *dVdt)
{
   
 for (int ii=0;ii<nCells;ii++) 
 {
   double Vm = VM[ii]; 
   double *states = States[ii].state; 
   int cellType = States[ii].cellType; 

   dVdt[ii]=0.0; 
   double itmp0,itmp5,itmp6 ; 
   {
     double fv1=fit[1]->afunc(Vm, fit[1]->aparms); 
     double fv2=fit[2]->afunc(Vm, fit[2]->aparms); 
     double x1 = states[Ca_i]*c26; 
     double x2 = SQ(states[Ca_i]*c27); 
     double x3 = SQ(states[Ca_i]*c28+c29); 
     double sigm1 = sigm(x1); 
     double sigm2 = sigm(x2); 
     double sigm3 = sigm(x3); 
     double dV3 = Vm- 0.5*c3*log(states[Ca_i]) - c8;

     itmp0 = CUBE(states[Na_i])*fv1-states[Ca_i]*fv2; 
     double itmp4 = c7*dV3+c24*sigm1; 
     itmp5=  c43*(states[Ca_i] -  states[Ca_SR])+c44*sigm2;      
     itmp6 = c23*(states[Ca_ss] - states[Ca_i]);
     states[Ca_i]    += dt*(sigm3*(itmp4-itmp0+itmp6*c15-itmp5*c16));
     dVdt[ii] -= itmp4*c10;
    }

   {

     double x0 = states[Na_i]*c25; 
     double sigm0 = sigm(x0); 

     double dV0 = Vm -c3*log(states[K_i]) -c5;
     double dV1 = Vm -c3*log(states[Na_i])-c4;
     double dV2 = dV0-c3*logSeries(c2*states[Na_i]/states[K_i])+c5-c6; // Assumption:  c2*states[Na_i]/states[K_i] is small; 

     double fv0=fit[0]->afunc(Vm, fit[0]->aparms); 
     double fv5=fit[5]->afunc(Vm, fit[5]->aparms); 
     double fv6=fit[6]->afunc(dV0,fit[6]->aparms); 
     double fd =  fv5  +  fv6; 

     double tmp0 =  fd +  cB[cellType]*states[r_gate]*states[s_gate]+ c11*states[Xr1_gate]*states[Xr2_gate] ;
     double tmp1 =  c20*CUBE(states[m_gate])*states[h_gate]*states[j_gate]+c21;
     double tmp2 =  cA[cellType]*SQ(states[Xs_gate]);
     double itmpA =  sigm0 * fv0;                          //Sigm0
     double itmp2 = itmp0 - 1.5*itmpA+tmp1*dV1; 
     double itmp3 = itmpA + tmp0*dV0 +tmp2*dV2; 

     states[K_i]     += dt*itmp3;
     states[Na_i]    += dt*(itmp2+2.0*itmp0);
     dVdt[ii]  += itmp2*c22+itmp3*c22;
   }

   {
     double fv3=fit[3]->afunc(Vm, fit[3]->aparms); 
     double fv4=fit[4]->afunc(Vm, fit[4]->aparms); 
     double x4 = SQ(states[Ca_SR]*c30); 
     double x5 = SQ(states[Ca_SR]*c31+c32); 
     double x6 = SQ(states[Ca_ss]*c33+c34); 
     double sigm4 = sigm(x4); 
     double sigm5 = sigm(x5); 
     double sigm6 = sigm(x6); 
     double tmp8  = c18+c19*sigm4; //Sigm4
     double tmp9  = tmp8*states[Ca_ss]+c36; 
     double itmp1 = states[d_gate]*states[f_gate]*states[f2_gate]*states[fCass_gate]*(fv4-states[Ca_ss]*fv3);
     double sigmCa_ss =   SQ(states[Ca_ss])/(tmp8*c17 + SQ(states[Ca_ss]));
     double itmp7 = sigmCa_ss*states[R_prime]*(states[Ca_SR] - states[Ca_ss]);

     double t1 = 1.0/(1.0+SQ(20*states[Ca_ss])); 
     double mhu = 0.600000*t1+0.4000000;
     double tauR =    1.0/(80.0*t1+2.0);

     states[Ca_ss]   += dt*(sigm6*(itmp6+itmp7*c14+itmp1*c13));    
     states[Ca_SR]   += dt*(sigm5*(itmp5-c40*itmp7));
     states[R_prime] += dt*(c36 - tmp9*states[R_prime]);
     states[fCass_gate] +=  dt*(mhu - states[fCass_gate])*tauR; 
     dVdt[ii] += itmp1;
   }
   }
}


// Update Gates; 
void computeGateRates(double dt, int nCells, const double *VM, TT06DevState* States)
{
   
 PADE **gatefit = fit + 7; 
 for (int ii=0;ii<nCells;ii++) 
 {
   double Vm = VM[ii]; 
   double *states = States[ii].state; 
   int cellType = States[ii].cellType; 
   double *gate = states+gateOffset; 
   
   double mhu,tauR; 

   int i=0;

   mhu=gatefit[2*i]  ->afunc(Vm,gatefit[2*i  ]->aparms);       //mGate
   tauR=gatefit[2*i+1]->afunc(Vm,gatefit[2*i+1]->aparms); 
   double grate=  (mhu - gate[i])*tauR; 
   grate  *= (1-exp(-dt*tauR))/(dt*tauR);                   //mtauR can be very large   ~1140.0 use Rush_Larson to integrate. 
   gate[i] += dt*grate; 
   i++;

   mhu=gatefit[2*i]  ->afunc(Vm,gatefit[2*i  ]->aparms);      // hGate 
   tauR=gatefit[2*i+1]->afunc(Vm,gatefit[2*i+1]->aparms); 
   gate[i] +=  dt*(mhu - gate[i])*tauR; 
   i++;

   //mhu=gatefit[2*i]  ->afunc(Vm,gatefit[2*i  ]->aparms);    //jGate        note hMhu and  jMhu are equal. 
   tauR=gatefit[2*i+1]->afunc(Vm,gatefit[2*i+1]->aparms); 
   gate[i] +=  dt*(mhu - gate[i])*tauR; 
   i++;

   for (int k=0;k<7;k++)                    // other Gates. 
   {
      mhu=gatefit[2*i]  ->afunc(Vm,gatefit[2*i  ]->aparms); 
      tauR=gatefit[2*i+1]->afunc(Vm,gatefit[2*i+1]->aparms); 
      gate[i] +=  dt*(mhu - gate[i])*tauR; 
      i++; 
   }
   int l = 2*(i+(cellType != 0)) ; 
   mhu=gatefit[ l ]->afunc(Vm,gatefit[l ]->aparms);  //Note  sMhu depends on celltype
   tauR=gatefit[l + 1 ]->afunc(Vm,gatefit[ l+1  ]->aparms);  //Note  sTauR depends on celltype
   gate[i] +=  dt*(mhu - gate[i])*tauR;     //sGate
}

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
   double t1 = 1400.00/sqrt(1.0+exp((5.0 - Vm)/6.0));
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
double hTauRMod(double Vm, void *parms) 
{ 
   double dtauR,ddtauR; 
   double tauR = TauRecipMod(Vm,hParms,&dtauR,&ddtauR); 
   return tauR;
}
double jTauR(double Vm, void *parms) 
{ 
   double t1  = (Vm < -40.0 ? (( ( - 25428.0*(exp(( 0.244400*Vm))) -  6.94800e-06*(exp(( - 0.0439100*Vm))))*(Vm+37.7800))/1.0)/(1.0+(exp(( 0.311000*(Vm+79.2300))))) : 0.0);
   double t2 = (Vm < -40.0 ? ( 0.0242400*(exp(( - 0.0105200*Vm))))/(1.0+(exp(( - 0.137800*(Vm+40.1400))))) : ( 0.600000*(exp(( 0.0570000*Vm))))/(1.0+(exp(( - 0.100000*(Vm+32.0))))));
   double tauR  = (t1+t2);
   return tauR;
}
double jTauRMod(double Vm, void *parms) 
{ 
   double dtauR,ddtauR; 
   double tauR = TauRecipMod(Vm,jParms,&dtauR,&ddtauR); 
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
 //double mhu = 0.670000/(1.0+(exp ((Vm+35.0)/7.0)) )+0.330000;
   double mhu = 0.670000/(1.0+(exp(((Vm+35.0)/7.0))))+0.330000;
   return mhu;
}
double f2TauR(double Vm, void *parms) 
{ 
   double tau =  562.0*exp(-SQ((Vm+27.0))/240.0)+31.0/(1.0+(exp(((25.0 - Vm)/10.0))))+80.0/(1.0+(exp(((Vm+30.0)/10.0))));
   double tauR = 1/tau; 
   return tauR;
}
double sMhu0(double Vm, void *parms) 
{ 
   double mhu = 1.0/(1.0+(exp(((Vm+28.0)/5.0))));
   return mhu;
}
double sMhu1(double Vm, void *parms) 
{ 
   double  mhu=1.00000/(1.00000+(exp(((Vm+20.0000)/5.00000))));
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
   double tau =  85.0000*(exp((- (pow((Vm+45.0000), 2.00000))/320.000)))+5.00000/(1.00000+(exp(((Vm - 20.0000)/5.00000))))+3.00000;
   double tauR = 1/tau; 
   return tauR;
}
