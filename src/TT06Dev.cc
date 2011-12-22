#include "TT06Dev.hh"
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include "TT06DevFit.hh"

using namespace std;
void initState(double *STATES,int cellType);

double TT06Dev::defaultState_[nStateVar];

TT06Dev::TT06Dev(int cellType)
{
   static bool initialized = false;
   cellType_ = cellType; 
   if (! initialized)
   {
      initialized = true;
      initState(defaultState_,cellType);
      initCnst();
      c9 = get_c9(); 
      //Approx(1301,32,0.005); 
   }
   for (unsigned ii=0; ii<nStateVar; ++ii) states_[ii] = defaultState_[ii];
}

/** returns dVm/dt for the reaction part only. */
double TT06Dev::calc(double dt, double Vm, double iStim, double states[nStateVar])
{
   double dVdt = computeUpdates(dt, Vm, states, cellType_);
   states[K_i] += iStim*c9*dt ;

   return dVdt;
}

/*
 * STATES[Vmembrane] is V in component membrane (millivolt).
 * STATES[K_i] is K_i in component potassium_dynamics (millimolar).
 * STATES[Na_i] is Na_i in component sodium_dynamics (millimolar).
 * STATES[Ca_i] is Ca_i in component calcium_dynamics (millimolar).
 * STATES[Xr1_gate] is Xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * STATES[Xr2_gate] is Xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * STATES[Xs_gate] is Xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * STATES[m_gate] is m in component fast_sodium_current_m_gate (dimensionless).
 * STATES[h_gate] is h in component fast_sodium_current_h_gate (dimensionless).
 * STATES[j_gate] is j in component fast_sodium_current_j_gate (dimensionless).
 * STATES[r_gate] is r in component transient_outward_current_r_gate (dimensionless).
 * STATES[d_gate] is d in component L_type_Ca_current_d_gate (dimensionless).
 * STATES[f_gate] is f in component L_type_Ca_current_f_gate (dimensionless).
 * STATES[f2_gate] is f2 in component L_type_Ca_current_f2_gate (dimensionless).
 * STATES[s_gate] is s in component transient_outward_current_s_gate (dimensionless).
 * STATES[fCass_gate] is fCass in component L_type_Ca_current_fCass_gate (dimensionless).
 * STATES[Ca_ss] is Ca_ss in component calcium_dynamics (millimolar).
 * STATES[Ca_SR] is Ca_SR in component calcium_dynamics (millimolar).
 * STATES[R_prime] is R_prime in component calcium_dynamics (dimensionless).
*/
