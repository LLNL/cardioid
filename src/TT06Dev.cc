#include "TT06Dev.hh"
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#define SQ(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

//enum { Vmembrane,K_i,Na_i,Ca_i,Xr1_gate,Xr2_gate,Xs_gate,m_gate,h_gate,j_gate,r_gate,d_gate,f_gate,f2_gate,fCass_gate,s_gate,Ca_ss,Ca_SR,R_prime}; 
using namespace std;

double TT06Dev::constants_[59];
double TT06Dev::defaultState_[19];
namespace
{
   void initState(double *STATES, int cellType);
   void initCnst();
   double computeUpdates(double dt, double* STATES,int cellType);
   void  update_Xr1Gate(double dt, double Vm, double *gate);
}
namespace
{
//   static double cnst[59]; 
   static double cA[3],cB[3]; 
   static double c1,c2,c3,c4,c5,c6,c7,c8,c9;
   static double c10,c11,c12,c13,c14,c15,c16,c17,c18,c19;
   static double c20,c21,c22,c23,c24,c25,c26,c27,c28,c29;
   static double c30,c31,c32,c33,c34,c36,c40,c43,c44;
   static double f1,f2,f3,f4,f5,f6,f7,f8,f9,f10; 
}
namespace 
{
void update_Xr1Gate(double dt, double Vm, double *gate)
{
   double mhu = 1.0/(1.0+(exp(((- 26.0 - Vm)/7.0))));
   double t1 = 450.0/(1.0+(exp(((- 45.0 - Vm)/10.0))));
   double t2 = 6.0/(1.0+(exp(((Vm+30.0)/11.5000))));
   double tau =  t1*t2;
   double rate = (mhu - *gate)/tau;
   *gate += rate*dt;
}
// Update Xr2 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
void update_Xr2Gate(double dt, double Vm, double *gate)
{
   double mhu = 1.0/(1.0+(exp(((Vm+88.0)/24.0))));
   double t1 = 3.0/(1.0+(exp(((- 60.0 - Vm)/20.0))));
   double t2 = 1.12000/(1.0+(exp(((Vm - 60.0)/20.0))));
   double tau =  t1*t2;
   double rate = (mhu - *gate)/tau;
   *gate += rate*dt;
}
void update_XsGate(double dt, double Vm, double *gate)
{
   double mhu = 1.0/(1.0+(exp(((- 5.0 - Vm)/14.0))));
   double t1 = 1400.00/ sqrt(1.0+exp((5.0 - Vm)/6.0));
   double t2 = 1.0/(1.0+(exp(((Vm - 35.0)/15.0))));
   double tau  =  t1*t2+80.0;
   double rate = (mhu - *gate)/tau;
   *gate += rate*dt;
}
void update_mGate(double dt, double Vm, double *gate)
{
   double mhu = 1.0/SQ(1.0+exp((-56.8600 - Vm)/9.03000));
   double t1  = 1.0/(1.0+(exp(((- 60.0 - Vm)/5.0))));
   double t2  = 0.100000/(1.0+(exp(((Vm+35.0)/5.0))))+0.100000/(1.0+(exp(((Vm - 50.0)/200.0))));
   double tau =  t1*t2;
   //double rate = (mhu - *gate)/tau;
   *gate = mhu - (mhu -*gate) *exp(-dt/tau);
}
void update_hGate(double dt, double Vm, double *gate)
{
   double mhu = 1.0/SQ((1.0+(exp(((Vm+71.5500)/7.43000)))));
   double t1  = (Vm<- 40.0 ?  0.0570000*(exp((- (Vm+80.0)/6.80000))) : 0.0);
   double t2  = (Vm<- 40.0 ?  2.70000*(exp(( 0.0790000*Vm)))+ 310000.*(exp(( 0.348500*Vm))) : 0.770000/( 0.130000*(1.0+(exp(((Vm+10.6600)/- 11.1000))))));
   double tau = 1.0/(t1+t2);
   double rate = (mhu - *gate)/tau;
   *gate += rate*dt;
}
void update_jGate(double dt, double Vm, double *gate)
{
   double mhu = 1.0/SQ((1.0+(exp(((Vm+71.5500)/7.43000)))));
   double t1  = (Vm<- 40.0 ? (( ( - 25428.0*(exp(( 0.244400*Vm))) -  6.94800e-06*(exp(( - 0.0439100*Vm))))*(Vm+37.7800))/1.0)/(1.0+(exp(( 0.311000*(Vm+79.2300))))) : 0.0);
   double t2 = (Vm<- 40.0 ? ( 0.0242400*(exp(( - 0.0105200*Vm))))/(1.0+(exp(( - 0.137800*(Vm+40.1400))))) : ( 0.600000*(exp(( 0.0570000*Vm))))/(1.0+(exp(( - 0.100000*(Vm+32.0))))));
   double tau  = 1.0/(t1+t2);
   double rate = (mhu - *gate)/tau;
   *gate += rate*dt;
}
void update_rGate(double dt, double Vm, double *gate)
{
   double mhu = 1.0/(1.0+(exp(((20.0 - Vm)/6.0))));
   double tau =  9.50000*(exp((- SQ((Vm+40.0)))/1800.00))+0.800000;
   double rate = (mhu - *gate)/tau;
   *gate += rate*dt;
}
void update_dGate(double dt, double Vm, double *gate)
{
   double mhu = 1.0/(1.0+(exp(((- 8.0 - Vm)/7.50000))));
   double t1  = 1.40000/(1.0+(exp(((- 35.0 - Vm)/13.0))))+0.250000;
   double t2 = 1.40000/(1.0+(exp(((Vm+5.0)/5.0))));
   double t3 = 1.0/(1.0+(exp(((50.0 - Vm)/20.0))));
   double tau =  t1*t2+t3;
   double rate = (mhu - *gate)/tau;
   *gate += rate*dt;
}
void update_fGate(double dt, double Vm, double *gate)
{
   double mhu = 1.0/(1.0+(exp(((Vm+20.0)/7.0))));
   double tau =  1102.50*(exp((- SQ(Vm+27.0)/225.0)))+200.0/(1.0+(exp(((13.0 - Vm)/10.0))))+180.0/(1.0+(exp(((Vm+30.0)/10.0))))+20.0;
   double rate = (mhu - *gate)/tau;
   *gate += rate*dt;
}
void update_f2Gate(double dt, double Vm, double *gate)
{
   double mhu = 0.670000/(1.0+(exp(((Vm+35.0)/7.0))))+0.330000;
   double tau =  562.0*exp(-SQ((Vm+27.0))/240.0)+31.0/(1.0+(exp(((25.0 - Vm)/10.0))))+80.0/(1.0+(exp(((Vm+30.0)/10.0))));
   double rate = (mhu - *gate)/tau;
   *gate += rate*dt;
}
void update_fCassGate(double dt, double Ca_SS, double *gate)
{
   double t1 = 1.0/(1.0+SQ(20*Ca_SS)); 
   double mhu = 0.600000*t1+0.4000000;
   double tau =    80.0*t1+2.0;
   double rate = (mhu - *gate)/tau;
   *gate += rate*dt;
}
void update_sGate(double dt, double Vm, double *gate, int cellType)
{

   if (cellType==0) 
   {
      double mhu = 1.0/(1.0+(exp(((Vm+28.0)/5.0))));
      double tau =  1000.0*(exp((- (pow((Vm+67.0), 2.0))/1000.0)))+8.0;
      double rate = (mhu - *gate)/tau;
      *gate += rate*dt;
   }
   else
   {
      double mhu = 1.0/(1.0+(exp(((Vm+20.0)/5.0))));
      double tau =  85.0*(exp((- (pow((Vm+45.0), 2.0))/320.0)))+5.0/(1.0+(exp(((Vm - 20.0)/5.0))))+3.0;
      double rate = (mhu - *gate)/tau;
      *gate += rate*dt;
   }
}
}

void fit(double dt)
{
double A[11],B[11]; 
FILE *file = fopen("file.input","w"); 
for (double Vm=-100;Vm<50;Vm+=0.02)
{
int i=0; 
double x; 
x = 0; update_Xr1Gate(dt, Vm, &x); A[i] = x; x = 1; update_Xr1Gate(dt, Vm, &x); B[i] = x-A[i]; i++; 
x = 0; update_Xr2Gate(dt, Vm, &x); A[i] = x; x = 1; update_Xr2Gate(dt, Vm, &x); B[i] = x-A[i]; i++; 
x = 0; update_XsGate(dt, Vm, &x);  A[i] = x; x = 1; update_XsGate(dt, Vm, &x); B[i] = x-A[i]; i++; 
x = 0; update_mGate(dt, Vm, &x);   A[i] = x; x = 1; update_mGate(dt, Vm, &x); B[i] = x-A[i]; i++; 
x = 0; update_hGate(dt, Vm, &x);   A[i] = x; x = 1; update_hGate(dt, Vm, &x); B[i] = x-A[i]; i++; 
x = 0; update_jGate(dt, Vm, &x);   A[i] = x; x = 1; update_jGate(dt, Vm, &x); B[i] = x-A[i]; i++; 
x = 0; update_rGate(dt, Vm, &x);   A[i] = x; x = 1; update_rGate(dt, Vm, &x); B[i] = x-A[i]; i++; 
x = 0; update_dGate(dt, Vm, &x);   A[i] = x; x = 1; update_dGate(dt, Vm, &x); B[i] = x-A[i]; i++; 
x = 0; update_fGate(dt, Vm, &x);   A[i] = x; x = 1; update_fGate(dt, Vm, &x); B[i] = x-A[i]; i++; 
x = 0; update_f2Gate(dt, Vm, &x);   A[i] = x; x = 1; update_f2Gate(dt, Vm, &x); B[i] = x-A[i]; i++; 
A[i]=B[i]=0.0; i++;
fprintf(file,"%f ", Vm); 
for (int j=0;j<i;j++) fprintf(file,"%e %e ", A[j],B[j]);  
fprintf(file,"\n"); 
}
/*
update_Xr2Gate(dt, Vm, STATES+5);
update_XsGate(dt, Vm, STATES+6);
update_mGate(dt, Vm, STATES+7);
update_hGate(dt, Vm, STATES+8);
update_jGate(dt, Vm, STATES+9);
update_rGate(dt, Vm, STATES+10);
update_dGate(dt, Vm, STATES+11);
update_fGate(dt, Vm, STATES+12);
update_f2Gate(dt, Vm, STATES+13);
update_fCassGate(dt, Ca_SS, STATES+14);
update_sGate(dt, Vm, STATES+15,cellType);
*/
}


   

TT06Dev::TT06Dev(int cellType)
{
   static bool initialized = false;
   cellType_ = cellType; 
   if (! initialized)
   {
      initialized = true;
      initState(defaultState_,cellType);
      initCnst();
      fit(0.02); 
   }
   for (unsigned ii=0; ii<19; ++ii) states_[ii] = defaultState_[ii];
}

/** returns dVm/dt for the reaction part only. */
double TT06Dev::calc(double dt, double Vm, double iStim, double states[19])
{
   states[0] = Vm;
   double dVdt = computeUpdates(dt, states, cellType_);
   states[1] += iStim*c9*dt ;

   return dVdt;
}


namespace
{
/*
   There are a total of 19 entries in state variable arrays.
   There are a total of 53 entries in the constant variable array.
 */
/*
 * STATES[0] is V in component membrane (millivolt).
 * STATES[1] is K_i in component potassium_dynamics (millimolar).
 * STATES[2] is Na_i in component sodium_dynamics (millimolar).
 * STATES[3] is Ca_i in component calcium_dynamics (millimolar).
 * STATES[4] is Xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * STATES[5] is Xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * STATES[6] is Xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * STATES[7] is m in component fast_sodium_current_m_gate (dimensionless).
 * STATES[8] is h in component fast_sodium_current_h_gate (dimensionless).
 * STATES[9] is j in component fast_sodium_current_j_gate (dimensionless).
 * STATES[10] is r in component transient_outward_current_r_gate (dimensionless).
 * STATES[11] is d in component L_type_Ca_current_d_gate (dimensionless).
 * STATES[12] is f in component L_type_Ca_current_f_gate (dimensionless).
 * STATES[13] is f2 in component L_type_Ca_current_f2_gate (dimensionless).
 * STATES[14] is fCass in component L_type_Ca_current_fCass_gate (dimensionless).
 * STATES[15] is s in component transient_outward_current_s_gate (dimensionless).
 * STATES[16] is Ca_ss in component calcium_dynamics (millimolar).
 * STATES[17] is Ca_SR in component calcium_dynamics (millimolar).
 * STATES[18] is R_prime in component calcium_dynamics (dimensionless).
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
void
initState(double *STATES,int cellType)
{
if (cellType == 0) 
{
STATES[0] = -86.709;
STATES[1] = 138.4;
STATES[2] = 10.355;
STATES[3] = 0.00013;
STATES[4] = 0.00448;
STATES[5] = 0.476;
STATES[6] = 0.0087;
STATES[7] = 0.00155;
STATES[8] = 0.7573;
STATES[9] = 0.7225;
STATES[10] = 2.235e-8;
STATES[11] = 3.164e-5;
STATES[12] = 0.8009;
STATES[13] = 0.9778;
STATES[14] = 0.9953;
STATES[15] = 0.3212;
STATES[16] = 0.00036;
STATES[17] = 3.715;
STATES[18] = 0.9068;
}

if (cellType == 1) 
{
STATES[0] = -85.423;
STATES[1] = 138.52;
STATES[2] = 10.132;
STATES[3] = 0.000153;
STATES[4] = 0.0165;
STATES[5] = 0.473;
STATES[6] = 0.0174;
STATES[7] = 0.00165;
STATES[8] = 0.749;
STATES[9] = 0.6788;
STATES[10] = 2.347e-8;
STATES[11] = 3.288e-5;
STATES[12] = 0.7026;
STATES[13] = 0.9526;
STATES[14] = 0.9942;
STATES[15] = 0.999998;
STATES[16] = 0.00042;
STATES[17] = 4.272;
STATES[18] = 0.8978;
}

if (cellType==2) 
{
STATES[0] = -85.23;
STATES[1] = 136.89;
STATES[2] = 8.604;
STATES[3] = 0.000126;
STATES[4] = 0.00621;
STATES[5] = 0.4712;
STATES[6] = 0.0095;
STATES[7] = 0.00172;
STATES[8] = 0.7444;
STATES[9] = 0.7045;
STATES[10] = 2.42e-8;
STATES[11] = 3.373e-5;
STATES[12] = 0.7888;
STATES[13] = 0.9755;
STATES[14] = 0.9953;
STATES[15] = 0.999998;
STATES[16] = 0.00036;
STATES[17] = 3.64;
STATES[18] = 0.9073;
}
}
void initCnst()
{
double cnst[59]; 
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

cnst[53] = 0.392;  //endo
cnst[54] = 0.098;  //mid
cnst[55] = 0.392;  //Epi

cnst[56] = 0.073;   //endo
cnst[57] = 0.294;   //mid
cnst[58] = 0.294;   //Epi

c1 = cnst[2]/(cnst[0]*cnst[1]); 
c2 = cnst[9]; 
c3 = -1/c1; 
c4 = -c3*log(cnst[11]);
c5 = -c3*log(cnst[10]);
c6 = -c3*log(cnst[10]+cnst[9]*cnst[11]);
//
c8 = -c3*log(cnst[12]);
c9 = -cnst[3]/(cnst[4]*cnst[2]);

cA[0] =  c9*cnst[53]; 
cA[1] =  c9*cnst[54]; 
cA[2] =  c9*cnst[55]; 

cB[0] =  c9*cnst[56]; 
cB[1] =  c9*cnst[57]; 
cB[2] =  c9*cnst[58]; 

c10= 1/(0.5*c9);
c7 =  (0.25*cnst[19]*c9);
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
f8 = 15.0*f7;
f9 = 4.0*cnst[12];
f10 = c9*cnst[32];
}
double computeUpdates(double dt, double* STATES, int cellType)
{


double Vm = STATES[0]; 
double Ca_SS = STATES[16]; 

double expV1 = exp(-f1*Vm); 
double expV2 = exp(f1*(30-2.0*Vm)); 
double expV5 = exp(f5*Vm); 
double expV = expV5*expV1; 
double fv1 = f2/(1.0+0.1245*exp(-0.1*Vm*f1)+0.0353*expV1);
double fv2 = 1.0/(f3+ f4*expV);
double fv3 = expV5*fv2 ; 
double fv4 = expV* fv2*f6 ; 
double fv5 = (f7*Vm-f8)/(1.0-expV2);
double fv6 = fv5 * expV2 *f9;
double fv7 = f10/(1.0+exp((25.0 - Vm)/5.98));

#define logSeries(x) ((x)*(1.0-0.5*(x)))
double dV0 = 1.0*Vm -c3*log(STATES[1]) -c5;
double dV1 = 1.0*Vm -c3*log(STATES[2]) -c4;
//double dV2 = 1.0*Vm -c3*log(STATES[1]+c2*STATES[2]) -c6;
double dV3 = 2.0*Vm- c3*log(STATES[3]) - c8;
//double ss = c2*(STATES[2]/STATES[1]);
//double dV2 =Vm-c3*log(STATES[1]) -c3*log(1+ss) -c6;
//double dV2 =Vm-c3*log(STATES[1])-c3*logSeries(c2*STATES[2]/STATES[1]) -c6;
double dV2 =dV0-c3*logSeries(c2*STATES[2]/STATES[1])+c5-c6;

double xx  =  (3.0*exp(0.0002*dV0 + 0.02)+exp(0.1*dV0 - 1.0))/(1.0+exp( -0.5*dV0))*(10.0+10*exp(0.06*dV0 -12.0));
double fdV0 =  c9*c12/(1.0+xx);

double fdV3 =  c7*dV3;


double x[8]; 
x[0] = STATES[2]*c25; 
x[1] = STATES[3]*c26; 
x[2] = SQ(STATES[3]*c27); 
x[3] = SQ(STATES[3]*c28+c29); 
x[4] = SQ(STATES[17]*c30); 
x[6] = SQ(STATES[17]*c31+c32); 
x[7] = SQ(STATES[16]*c33+c34); 
#define sigm(x)   ((x)/(1.0+(x)) )
double sigm0 = sigm(x[0]); 
double sigm1 = sigm(x[1]); 
double sigm2 = sigm(x[2]); 
double sigm3 = sigm(x[3]); 
double sigm4 = sigm(x[4]); 
double sigm6 = sigm(x[6]); 
double sigm7 = sigm(x[7]); 


double tmp5 =    fdV0 +      cB[cellType]*STATES[10]*STATES[15]+ c11*STATES[4]*STATES[5] +  fv7;
double tmp6 =  c20*CUBE(STATES[7])*STATES[8]*STATES[9]+c21;
double tmp7 =  cA[cellType]*SQ(STATES[6]);
double tmp8  = c18+c19*sigm4; //Sigm4

double sigm5 =   SQ(STATES[16])/(tmp8*c17 + SQ(STATES[16]));
double  tmp9 =    tmp8*STATES[16]+c36; 

double itmpA =  sigm0 * fv1;                          //Sigm0
double itmp0 = CUBE(STATES[2])*fv3-STATES[3]*fv4 ; 
double itmp1 =  itmp0-1.5*itmpA+tmp6*dV1; 
double itmp2 = itmpA + tmp5*dV0+tmp7*dV2; 
double itmp3 = STATES[11]*STATES[12]*STATES[13]*STATES[14]*(fv6-STATES[16]*fv5);
double itmp4 = fdV3+c24*sigm1; 
double itmp5=  c43*(STATES[3] - STATES[17])+c44*sigm2;      
double itmp6 = c23*(STATES[16] - STATES[3]);
double itmp7 = sigm5*STATES[18]*(STATES[17] - STATES[16]);

double dVdt  = itmp3+itmp1*c22+itmp2*c22-itmp4*c10;
STATES[1]  += dt*(itmp2);
STATES[2]  += dt*(itmp1+2.0*itmp0);
STATES[3]  += dt*(sigm3*(itmp4-itmp0+itmp6*c15-itmp5*c16));
STATES[16] += dt*(sigm7*(itmp6+itmp7*c14+itmp3*c13));    
STATES[17] += dt*(sigm6*(itmp5-c40*itmp7));
STATES[18] += dt*(c36 - tmp9*STATES[18]);

update_Xr1Gate(dt, Vm, STATES+4);
update_Xr2Gate(dt, Vm, STATES+5);
update_XsGate(dt, Vm, STATES+6);
update_mGate(dt, Vm, STATES+7);
update_hGate(dt, Vm, STATES+8);
update_jGate(dt, Vm, STATES+9);
update_rGate(dt, Vm, STATES+10);
update_dGate(dt, Vm, STATES+11);
update_fGate(dt, Vm, STATES+12);
update_f2Gate(dt, Vm, STATES+13);
update_fCassGate(dt, Ca_SS, STATES+14);
update_sGate(dt, Vm, STATES+15,cellType);

/*
static FILE *file = NULL; 
static double time = 0.0; 
static int loop =0; 
static double xmin[7],xmax[7]; 
if (file == NULL) 
{	file = fopen("info.data","w"); 
//	for (int i=0;i<7;i++)  xmin[i]=xmax[i] = x[i]; 
}
if (loop % 10 ==0) fprintf(file,"%14.10f %e %24.14e %24.14e %24.14e %24.14e %24.14e %24.14e\n",time,RATES[0],RATES[1],RATES[2],RATES[3],RATES[16],RATES[17],STATES[18]); 
if (loop % 1000 ==0) fprintf(file,"    %f %e %e %e %e %e %e %e\n",time,x[0],x[1],x[2],x[3],x[4],x[5],xmin[6]); 
for (int i=0;i<7;i++) 
{
	if (x[i] < xmin[i] ) xmin[i] = x[i]; 
	if (x[i] > xmax[i] ) xmax[i] = x[i]; 
}
if (loop % 1000 ==0) fprintf(file,"min %f %e %e %e %e %e %e %e\n",time,xmin[0],xmin[1],xmin[2],xmin[3],xmin[4],xmin[5],xmin[6]); 
if (loop % 1000 ==0) fprintf(file,"max %f %e %e %e %e %e %e %e\n",time,xmax[0],xmax[1],xmax[2],xmax[3],xmax[4],xmax[5],xmax[6]); 
time += 0.02; 
loop++; 
*/

return dVdt; 
}
}
