#ifndef TT06FUNC_HH
#define TT06FUNC_HH
#include "pade.hh" 
#define SQ(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

/*
* RATES[0] is d/dt V in component membrane (millivolt).
 * RATES[4] is d/dt Xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * RATES[5] is d/dt Xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * RATES[6] is d/dt Xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * RATES[7] is d/dt m in component fast_sodium_current_m_gate (dimensionless).
 * RATES[8] is d/dt h in component fast_sodium_current_h_gate (dimensionless).
 * RATES[9] is d/dt j in component fast_sodium_current_j_gate (dimensionless).
 * RATES[11] is d/dt d in component L_type_Ca_current_d_gate (dimensionless).
 * RATES[12] is d/dt f in component L_type_Ca_current_f_gate (dimensionless).
 * RATES[13] is d/dt f2 in component L_type_Ca_current_f2_gate (dimensionless).
 * RATES[14] is d/dt fCass in component L_type_Ca_current_fCass_gate (dimensionless).
 * RATES[15] is d/dt s in component transient_outward_current_s_gate (dimensionless).
 * RATES[16] is d/dt r in component transient_outward_current_r_gate (dimensionless).
 * RATES[18] is d/dt R_prime in component calcium_dynamics (dimensionless).
 * RATES[3] is d/dt Ca_i in component calcium_dynamics (millimolar).
 * RATES[17] is d/dt Ca_SR in component calcium_dynamics (millimolar).
 * RATES[10] is d/dt Ca_ss in component calcium_dynamics (millimolar).
 * RATES[2] is d/dt Na_i in component sodium_dynamics (millimolar).
 * RATES[1] is d/dt K_i in component potassium_dynamics (millimolar).
*/
// int mapCell2Dev[]                  {1,2,3,10,17,18,14,7,8,9,4,5,6,16,11,12,13,15};
namespace TT06Func
{
//enum TT06STATE { Ca_i, dVK_i, Na_i, Ca_ss, Ca_SR, R_prime, fCass_gate, m_gate, h_gate, j_gate, Xr1_gate, Xr2_gate, Xs_gate, r_gate, d_gate, f_gate, f2_gate,  jL_gate, s_gate, nStateVar} ; 
enum TT06STATE { Ca_i, dVK_i, Na_i, Ca_ss, Ca_SR, R_prime, fCass_gate, nStateVar} ; 
enum TT06GATES { m_gateN, h_gateN, j_gateN, Xr1_gateN, Xr2_gateN, Xs_gateN, r_gateN, d_gateN, f_gateN, f2_gateN,  jL_gateN, s_gateN, nGatesVar} ; 

struct TT06DevState
{
   int cellType; 
   double P_NaK,g_Ks,g_to,g_NaL; 
   double state[nStateVar];
};

#define gateOffset 7
void initState(TT06DevState *cell,double *gate, int cellType);
void initGate(double **gates, int ii, double *gateI);
void initCnst();
void updateNonGate(double dt, int n, const double *Vm, TT06DevState *cell, int offset, double **gates, double *dVdt);
void updateGate(double dt, int n, const double *Vm, TT06DevState *cell, int offset, double **gates);
void updateGate0(double dt, int n, const double *Vm, TT06DevState *cell, int offset, double **gates);
void updateGate1(double dt, int n, const double *Vm, TT06DevState *cell, int offset, double **gates);
void updateGate2(double dt, int n, const double *Vm, TT06DevState *cell, int offset, double **gates);
void updateGate3(double dt, int n, const double *Vm, TT06DevState *cell, int offset, double **gates);
void updateGateOLD(double dt, int n, const double *Vm, TT06DevState *cell);
double get_c9();
PADE **makeFit(double tol, double V0, double V1, double deltaV, int mod);
void writeFit(PADE **fit); 
double defaultVoltage(int cellType);
};

#endif
