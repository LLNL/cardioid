#include "TT06_CellML_Mid.hh"
#include <cmath>

using namespace std;

double TT06_CellML_Mid::constants_[53];
double TT06_CellML_Mid::defaultState_[19];

namespace
{
   void
      initConsts(double* CONSTANTS, double* RATES, double *STATES);
   void
      computeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC);
}

TT06_CellML_Mid::TT06_CellML_Mid()
{
   static bool initialized = false;
   if (! initialized)
   {
      initialized = true;
      double dummy[19];
      initConsts(constants_, dummy, defaultState_);
   }
   for (unsigned ii=0; ii<19; ++ii)
      states_[ii] = defaultState_[ii];
}

/** returns dVm/dt for the reaction part only. */
double TT06_CellML_Mid::calc(double dt, double Vm, double iStim)
{
   double dummy  = 0.;
   double rates[19];
   double algebraic[70];
   algebraic[12] = iStim;
   states_[0] = Vm;
   computeRates(dummy, constants_, rates, states_, algebraic);

   for (unsigned ii=1; ii<19; ++ii)
      states_[ii] += rates[ii] * dt;
   
   return rates[0]+iStim;
}


/** Everything below this point is the TT06 mid model generated code
 *  from the CellML web site (retrieved 8-Nov-2011).  This code has been
 *  altered as follows:
 *
 *  1.  Placed in the anonymous namespace so that it will not conflict
 *  with other compilation units.
 *
 *  2. The computeVariables function was deleted since it is not needed
 *  for a forward Euler method
 *
 *  3. The computation of ALGEBRAIC[12] in computeRates is commented
 *  out.  ALGEBRAIC[12] corresponds to the stimulus current.  It is
 *  expected that ALGEBRAIC[12] will be loaded with a physically relevant
 *  value by the caller.  Note that VOI is used only in the calculation
 *  of ALGEBRAIC[12] so it effectively becomes an unused parameter as do
 *  CONSTANTS[5:8] that specify the waveform of the stimulus.
 *
 *  Note that we have not modified the calculation of RATE[0].  This
 *  means it still includes the contribution due to iStim (aka
 *  ALGEBRAIC[12]).  This is not the convention that we use in our code.
 *  We expect that the reaction model will return only the dVm/dt due to
 *  reaction part only.  The external current contribution will be added
 *  and integrated in the main loop.  We have not changed the code so as
 *  to use the lightest possible touch on the CellML code. */

namespace
{
/*
   There are a total of 70 entries in the algebraic variable array.
   There are a total of 19 entries in each of the rate and state variable arrays.
   There are a total of 53 entries in the constant variable array.
 */
/*
 * VOI is time in component environment (millisecond).
 * STATES[0] is V in component membrane (millivolt).
 * CONSTANTS[0] is R in component membrane (joule_per_mole_kelvin).
 * CONSTANTS[1] is T in component membrane (kelvin).
 * CONSTANTS[2] is F in component membrane (coulomb_per_millimole).
 * CONSTANTS[3] is Cm in component membrane (microF).
 * CONSTANTS[4] is V_c in component membrane (micrometre3).
 * ALGEBRAIC[47] is i_K1 in component inward_rectifier_potassium_current (picoA_per_picoF).
 * ALGEBRAIC[54] is i_to in component transient_outward_current (picoA_per_picoF).
 * ALGEBRAIC[48] is i_Kr in component rapid_time_dependent_potassium_current (picoA_per_picoF).
 * ALGEBRAIC[49] is i_Ks in component slow_time_dependent_potassium_current (picoA_per_picoF).
 * ALGEBRAIC[52] is i_CaL in component L_type_Ca_current (picoA_per_picoF).
 * ALGEBRAIC[55] is i_NaK in component sodium_potassium_pump_current (picoA_per_picoF).
 * ALGEBRAIC[50] is i_Na in component fast_sodium_current (picoA_per_picoF).
 * ALGEBRAIC[51] is i_b_Na in component sodium_background_current (picoA_per_picoF).
 * ALGEBRAIC[56] is i_NaCa in component sodium_calcium_exchanger_current (picoA_per_picoF).
 * ALGEBRAIC[53] is i_b_Ca in component calcium_background_current (picoA_per_picoF).
 * ALGEBRAIC[58] is i_p_K in component potassium_pump_current (picoA_per_picoF).
 * ALGEBRAIC[57] is i_p_Ca in component calcium_pump_current (picoA_per_picoF).
 * ALGEBRAIC[12] is i_Stim in component membrane (picoA_per_picoF).
 * CONSTANTS[5] is stim_start in component membrane (millisecond).
 * CONSTANTS[6] is stim_period in component membrane (millisecond).
 * CONSTANTS[7] is stim_duration in component membrane (millisecond).
 * CONSTANTS[8] is stim_amplitude in component membrane (picoA_per_picoF).
 * ALGEBRAIC[25] is E_Na in component reversal_potentials (millivolt).
 * ALGEBRAIC[33] is E_K in component reversal_potentials (millivolt).
 * ALGEBRAIC[41] is E_Ks in component reversal_potentials (millivolt).
 * ALGEBRAIC[43] is E_Ca in component reversal_potentials (millivolt).
 * CONSTANTS[9] is P_kna in component reversal_potentials (dimensionless).
 * CONSTANTS[10] is K_o in component potassium_dynamics (millimolar).
 * CONSTANTS[11] is Na_o in component sodium_dynamics (millimolar).
 * STATES[1] is K_i in component potassium_dynamics (millimolar).
 * STATES[2] is Na_i in component sodium_dynamics (millimolar).
 * CONSTANTS[12] is Ca_o in component calcium_dynamics (millimolar).
 * STATES[3] is Ca_i in component calcium_dynamics (millimolar).
 * CONSTANTS[13] is g_K1 in component inward_rectifier_potassium_current (nanoS_per_picoF).
 * ALGEBRAIC[46] is xK1_inf in component inward_rectifier_potassium_current (dimensionless).
 * ALGEBRAIC[44] is alpha_K1 in component inward_rectifier_potassium_current (dimensionless).
 * ALGEBRAIC[45] is beta_K1 in component inward_rectifier_potassium_current (dimensionless).
 * CONSTANTS[14] is g_Kr in component rapid_time_dependent_potassium_current (nanoS_per_picoF).
 * STATES[4] is Xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * STATES[5] is Xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * ALGEBRAIC[0] is xr1_inf in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * ALGEBRAIC[13] is alpha_xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * ALGEBRAIC[26] is beta_xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * ALGEBRAIC[34] is tau_xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (millisecond).
 * ALGEBRAIC[1] is xr2_inf in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * ALGEBRAIC[14] is alpha_xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * ALGEBRAIC[27] is beta_xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * ALGEBRAIC[35] is tau_xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (millisecond).
 * CONSTANTS[15] is g_Ks in component slow_time_dependent_potassium_current (nanoS_per_picoF).
 * STATES[6] is Xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * ALGEBRAIC[2] is xs_inf in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * ALGEBRAIC[15] is alpha_xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * ALGEBRAIC[28] is beta_xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * ALGEBRAIC[36] is tau_xs in component slow_time_dependent_potassium_current_Xs_gate (millisecond).
 * CONSTANTS[16] is g_Na in component fast_sodium_current (nanoS_per_picoF).
 * STATES[7] is m in component fast_sodium_current_m_gate (dimensionless).
 * STATES[8] is h in component fast_sodium_current_h_gate (dimensionless).
 * STATES[9] is j in component fast_sodium_current_j_gate (dimensionless).
 * ALGEBRAIC[3] is m_inf in component fast_sodium_current_m_gate (dimensionless).
 * ALGEBRAIC[16] is alpha_m in component fast_sodium_current_m_gate (dimensionless).
 * ALGEBRAIC[29] is beta_m in component fast_sodium_current_m_gate (dimensionless).
 * ALGEBRAIC[37] is tau_m in component fast_sodium_current_m_gate (millisecond).
 * ALGEBRAIC[4] is h_inf in component fast_sodium_current_h_gate (dimensionless).
 * ALGEBRAIC[17] is alpha_h in component fast_sodium_current_h_gate (per_millisecond).
 * ALGEBRAIC[30] is beta_h in component fast_sodium_current_h_gate (per_millisecond).
 * ALGEBRAIC[38] is tau_h in component fast_sodium_current_h_gate (millisecond).
 * ALGEBRAIC[5] is j_inf in component fast_sodium_current_j_gate (dimensionless).
 * ALGEBRAIC[18] is alpha_j in component fast_sodium_current_j_gate (per_millisecond).
 * ALGEBRAIC[31] is beta_j in component fast_sodium_current_j_gate (per_millisecond).
 * ALGEBRAIC[39] is tau_j in component fast_sodium_current_j_gate (millisecond).
 * CONSTANTS[17] is g_bna in component sodium_background_current (nanoS_per_picoF).
 * CONSTANTS[18] is g_CaL in component L_type_Ca_current (litre_per_farad_second).
 * STATES[10] is Ca_ss in component calcium_dynamics (millimolar).
 * STATES[11] is d in component L_type_Ca_current_d_gate (dimensionless).
 * STATES[12] is f in component L_type_Ca_current_f_gate (dimensionless).
 * STATES[13] is f2 in component L_type_Ca_current_f2_gate (dimensionless).
 * STATES[14] is fCass in component L_type_Ca_current_fCass_gate (dimensionless).
 * ALGEBRAIC[6] is d_inf in component L_type_Ca_current_d_gate (dimensionless).
 * ALGEBRAIC[19] is alpha_d in component L_type_Ca_current_d_gate (dimensionless).
 * ALGEBRAIC[32] is beta_d in component L_type_Ca_current_d_gate (dimensionless).
 * ALGEBRAIC[40] is gamma_d in component L_type_Ca_current_d_gate (millisecond).
 * ALGEBRAIC[42] is tau_d in component L_type_Ca_current_d_gate (millisecond).
 * ALGEBRAIC[7] is f_inf in component L_type_Ca_current_f_gate (dimensionless).
 * ALGEBRAIC[20] is tau_f in component L_type_Ca_current_f_gate (millisecond).
 * ALGEBRAIC[8] is f2_inf in component L_type_Ca_current_f2_gate (dimensionless).
 * ALGEBRAIC[21] is tau_f2 in component L_type_Ca_current_f2_gate (millisecond).
 * ALGEBRAIC[9] is fCass_inf in component L_type_Ca_current_fCass_gate (dimensionless).
 * ALGEBRAIC[22] is tau_fCass in component L_type_Ca_current_fCass_gate (millisecond).
 * CONSTANTS[19] is g_bca in component calcium_background_current (nanoS_per_picoF).
 * CONSTANTS[20] is g_to in component transient_outward_current (nanoS_per_picoF).
 * STATES[15] is s in component transient_outward_current_s_gate (dimensionless).
 * STATES[16] is r in component transient_outward_current_r_gate (dimensionless).
 * ALGEBRAIC[10] is s_inf in component transient_outward_current_s_gate (dimensionless).
 * ALGEBRAIC[23] is tau_s in component transient_outward_current_s_gate (millisecond).
 * ALGEBRAIC[11] is r_inf in component transient_outward_current_r_gate (dimensionless).
 * ALGEBRAIC[24] is tau_r in component transient_outward_current_r_gate (millisecond).
 * CONSTANTS[21] is P_NaK in component sodium_potassium_pump_current (picoA_per_picoF).
 * CONSTANTS[22] is K_mk in component sodium_potassium_pump_current (millimolar).
 * CONSTANTS[23] is K_mNa in component sodium_potassium_pump_current (millimolar).
 * CONSTANTS[24] is K_NaCa in component sodium_calcium_exchanger_current (picoA_per_picoF).
 * CONSTANTS[25] is K_sat in component sodium_calcium_exchanger_current (dimensionless).
 * CONSTANTS[26] is alpha in component sodium_calcium_exchanger_current (dimensionless).
 * CONSTANTS[27] is gamma in component sodium_calcium_exchanger_current (dimensionless).
 * CONSTANTS[28] is Km_Ca in component sodium_calcium_exchanger_current (millimolar).
 * CONSTANTS[29] is Km_Nai in component sodium_calcium_exchanger_current (millimolar).
 * CONSTANTS[30] is g_pCa in component calcium_pump_current (picoA_per_picoF).
 * CONSTANTS[31] is K_pCa in component calcium_pump_current (millimolar).
 * CONSTANTS[32] is g_pK in component potassium_pump_current (nanoS_per_picoF).
 * STATES[17] is Ca_SR in component calcium_dynamics (millimolar).
 * ALGEBRAIC[67] is i_rel in component calcium_dynamics (millimolar_per_millisecond).
 * ALGEBRAIC[59] is i_up in component calcium_dynamics (millimolar_per_millisecond).
 * ALGEBRAIC[60] is i_leak in component calcium_dynamics (millimolar_per_millisecond).
 * ALGEBRAIC[61] is i_xfer in component calcium_dynamics (millimolar_per_millisecond).
 * ALGEBRAIC[66] is O in component calcium_dynamics (dimensionless).
 * STATES[18] is R_prime in component calcium_dynamics (dimensionless).
 * ALGEBRAIC[64] is k1 in component calcium_dynamics (per_millimolar2_per_millisecond).
 * ALGEBRAIC[65] is k2 in component calcium_dynamics (per_millimolar_per_millisecond).
 * CONSTANTS[33] is k1_prime in component calcium_dynamics (per_millimolar2_per_millisecond).
 * CONSTANTS[34] is k2_prime in component calcium_dynamics (per_millimolar_per_millisecond).
 * CONSTANTS[35] is k3 in component calcium_dynamics (per_millisecond).
 * CONSTANTS[36] is k4 in component calcium_dynamics (per_millisecond).
 * CONSTANTS[37] is EC in component calcium_dynamics (millimolar).
 * CONSTANTS[38] is max_sr in component calcium_dynamics (dimensionless).
 * CONSTANTS[39] is min_sr in component calcium_dynamics (dimensionless).
 * ALGEBRAIC[62] is kcasr in component calcium_dynamics (dimensionless).
 * CONSTANTS[40] is V_rel in component calcium_dynamics (per_millisecond).
 * CONSTANTS[41] is V_xfer in component calcium_dynamics (per_millisecond).
 * CONSTANTS[42] is K_up in component calcium_dynamics (millimolar).
 * CONSTANTS[43] is V_leak in component calcium_dynamics (per_millisecond).
 * CONSTANTS[44] is Vmax_up in component calcium_dynamics (millimolar_per_millisecond).
 * ALGEBRAIC[63] is Ca_i_bufc in component calcium_dynamics (dimensionless).
 * ALGEBRAIC[68] is Ca_sr_bufsr in component calcium_dynamics (dimensionless).
 * ALGEBRAIC[69] is Ca_ss_bufss in component calcium_dynamics (dimensionless).
 * CONSTANTS[45] is Buf_c in component calcium_dynamics (millimolar).
 * CONSTANTS[46] is K_buf_c in component calcium_dynamics (millimolar).
 * CONSTANTS[47] is Buf_sr in component calcium_dynamics (millimolar).
 * CONSTANTS[48] is K_buf_sr in component calcium_dynamics (millimolar).
 * CONSTANTS[49] is Buf_ss in component calcium_dynamics (millimolar).
 * CONSTANTS[50] is K_buf_ss in component calcium_dynamics (millimolar).
 * CONSTANTS[51] is V_sr in component calcium_dynamics (micrometre3).
 * CONSTANTS[52] is V_ss in component calcium_dynamics (micrometre3).
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
void
initConsts(double* CONSTANTS, double* RATES, double *STATES)
{
STATES[0] = -85.423;
CONSTANTS[0] = 8314.472;
CONSTANTS[1] = 310;
CONSTANTS[2] = 96485.3415;
CONSTANTS[3] = 0.185;
CONSTANTS[4] = 0.016404;
CONSTANTS[5] = 10;
CONSTANTS[6] = 1000;
CONSTANTS[7] = 1;
CONSTANTS[8] = 52;
CONSTANTS[9] = 0.03;
CONSTANTS[10] = 5.4;
CONSTANTS[11] = 140;
STATES[1] = 138.52;
STATES[2] = 10.132;
CONSTANTS[12] = 2;
STATES[3] = 0.000153;
CONSTANTS[13] = 5.405;
CONSTANTS[14] = 0.153;
STATES[4] = 0.0165;
STATES[5] = 0.473;
CONSTANTS[15] = 0.098;
STATES[6] = 0.0174;
CONSTANTS[16] = 14.838;
STATES[7] = 0.00165;
STATES[8] = 0.749;
STATES[9] = 0.6788;
CONSTANTS[17] = 0.00029;
CONSTANTS[18] = 0.0000398;
STATES[10] = 0.00042;
STATES[11] = 3.288e-5;
STATES[12] = 0.7026;
STATES[13] = 0.9526;
STATES[14] = 0.9942;
CONSTANTS[19] = 0.000592;
CONSTANTS[20] = 0.294;
STATES[15] = 0.999998;
STATES[16] = 2.347e-8;
CONSTANTS[21] = 2.724;
CONSTANTS[22] = 1;
CONSTANTS[23] = 40;
CONSTANTS[24] = 1000;
CONSTANTS[25] = 0.1;
CONSTANTS[26] = 2.5;
CONSTANTS[27] = 0.35;
CONSTANTS[28] = 1.38;
CONSTANTS[29] = 87.5;
CONSTANTS[30] = 0.1238;
CONSTANTS[31] = 0.0005;
CONSTANTS[32] = 0.0146;
STATES[17] = 4.272;
STATES[18] = 0.8978;
CONSTANTS[33] = 0.15;
CONSTANTS[34] = 0.045;
CONSTANTS[35] = 0.06;
CONSTANTS[36] = 0.005;
CONSTANTS[37] = 1.5;
CONSTANTS[38] = 2.5;
CONSTANTS[39] = 1;
CONSTANTS[40] = 0.102;
CONSTANTS[41] = 0.0038;
CONSTANTS[42] = 0.00025;
CONSTANTS[43] = 0.00036;
CONSTANTS[44] = 0.006375;
CONSTANTS[45] = 0.2;
CONSTANTS[46] = 0.001;
CONSTANTS[47] = 10;
CONSTANTS[48] = 0.3;
CONSTANTS[49] = 0.4;
CONSTANTS[50] = 0.00025;
CONSTANTS[51] = 0.001094;
CONSTANTS[52] = 0.00005468;
}
void
computeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
ALGEBRAIC[7] = 1.00000/(1.00000+(exp(((STATES[0]+20.0000)/7.00000))));
ALGEBRAIC[20] =  1102.50*(exp((- (pow((STATES[0]+27.0000), 2.00000))/225.000)))+200.000/(1.00000+(exp(((13.0000 - STATES[0])/10.0000))))+180.000/(1.00000+(exp(((STATES[0]+30.0000)/10.0000))))+20.0000;
RATES[12] = (ALGEBRAIC[7] - STATES[12])/ALGEBRAIC[20];
ALGEBRAIC[8] = 0.670000/(1.00000+(exp(((STATES[0]+35.0000)/7.00000))))+0.330000;
ALGEBRAIC[21] =  562.000*(exp((- (pow((STATES[0]+27.0000), 2.00000))/240.000)))+31.0000/(1.00000+(exp(((25.0000 - STATES[0])/10.0000))))+80.0000/(1.00000+(exp(((STATES[0]+30.0000)/10.0000))));
RATES[13] = (ALGEBRAIC[8] - STATES[13])/ALGEBRAIC[21];
ALGEBRAIC[9] = 0.600000/(1.00000+(pow((STATES[10]/0.0500000), 2.00000)))+0.400000;
ALGEBRAIC[22] = 80.0000/(1.00000+(pow((STATES[10]/0.0500000), 2.00000)))+2.00000;
RATES[14] = (ALGEBRAIC[9] - STATES[14])/ALGEBRAIC[22];
ALGEBRAIC[10] = 1.00000/(1.00000+(exp(((STATES[0]+20.0000)/5.00000))));
ALGEBRAIC[23] =  85.0000*(exp((- (pow((STATES[0]+45.0000), 2.00000))/320.000)))+5.00000/(1.00000+(exp(((STATES[0] - 20.0000)/5.00000))))+3.00000;
RATES[15] = (ALGEBRAIC[10] - STATES[15])/ALGEBRAIC[23];
ALGEBRAIC[11] = 1.00000/(1.00000+(exp(((20.0000 - STATES[0])/6.00000))));
ALGEBRAIC[24] =  9.50000*(exp((- (pow((STATES[0]+40.0000), 2.00000))/1800.00)))+0.800000;
RATES[16] = (ALGEBRAIC[11] - STATES[16])/ALGEBRAIC[24];
ALGEBRAIC[0] = 1.00000/(1.00000+(exp(((- 26.0000 - STATES[0])/7.00000))));
ALGEBRAIC[13] = 450.000/(1.00000+(exp(((- 45.0000 - STATES[0])/10.0000))));
ALGEBRAIC[26] = 6.00000/(1.00000+(exp(((STATES[0]+30.0000)/11.5000))));
ALGEBRAIC[34] =  1.00000*ALGEBRAIC[13]*ALGEBRAIC[26];
RATES[4] = (ALGEBRAIC[0] - STATES[4])/ALGEBRAIC[34];
ALGEBRAIC[1] = 1.00000/(1.00000+(exp(((STATES[0]+88.0000)/24.0000))));
ALGEBRAIC[14] = 3.00000/(1.00000+(exp(((- 60.0000 - STATES[0])/20.0000))));
ALGEBRAIC[27] = 1.12000/(1.00000+(exp(((STATES[0] - 60.0000)/20.0000))));
ALGEBRAIC[35] =  1.00000*ALGEBRAIC[14]*ALGEBRAIC[27];
RATES[5] = (ALGEBRAIC[1] - STATES[5])/ALGEBRAIC[35];
ALGEBRAIC[2] = 1.00000/(1.00000+(exp(((- 5.00000 - STATES[0])/14.0000))));
ALGEBRAIC[15] = 1400.00/ pow((1.00000+(exp(((5.00000 - STATES[0])/6.00000)))), 1.0 / 2);
ALGEBRAIC[28] = 1.00000/(1.00000+(exp(((STATES[0] - 35.0000)/15.0000))));
ALGEBRAIC[36] =  1.00000*ALGEBRAIC[15]*ALGEBRAIC[28]+80.0000;
RATES[6] = (ALGEBRAIC[2] - STATES[6])/ALGEBRAIC[36];
ALGEBRAIC[3] = 1.00000/(pow((1.00000+(exp(((- 56.8600 - STATES[0])/9.03000)))), 2.00000));
ALGEBRAIC[16] = 1.00000/(1.00000+(exp(((- 60.0000 - STATES[0])/5.00000))));
ALGEBRAIC[29] = 0.100000/(1.00000+(exp(((STATES[0]+35.0000)/5.00000))))+0.100000/(1.00000+(exp(((STATES[0] - 50.0000)/200.000))));
ALGEBRAIC[37] =  1.00000*ALGEBRAIC[16]*ALGEBRAIC[29];
RATES[7] = (ALGEBRAIC[3] - STATES[7])/ALGEBRAIC[37];
ALGEBRAIC[4] = 1.00000/(pow((1.00000+(exp(((STATES[0]+71.5500)/7.43000)))), 2.00000));
ALGEBRAIC[17] = (STATES[0]<- 40.0000 ?  0.0570000*(exp((- (STATES[0]+80.0000)/6.80000))) : 0.00000);
ALGEBRAIC[30] = (STATES[0]<- 40.0000 ?  2.70000*(exp(( 0.0790000*STATES[0])))+ 310000.*(exp(( 0.348500*STATES[0]))) : 0.770000/( 0.130000*(1.00000+(exp(((STATES[0]+10.6600)/- 11.1000))))));
ALGEBRAIC[38] = 1.00000/(ALGEBRAIC[17]+ALGEBRAIC[30]);
RATES[8] = (ALGEBRAIC[4] - STATES[8])/ALGEBRAIC[38];
ALGEBRAIC[5] = 1.00000/(pow((1.00000+(exp(((STATES[0]+71.5500)/7.43000)))), 2.00000));
ALGEBRAIC[18] = (STATES[0]<- 40.0000 ? (( ( - 25428.0*(exp(( 0.244400*STATES[0]))) -  6.94800e-06*(exp(( - 0.0439100*STATES[0]))))*(STATES[0]+37.7800))/1.00000)/(1.00000+(exp(( 0.311000*(STATES[0]+79.2300))))) : 0.00000);
ALGEBRAIC[31] = (STATES[0]<- 40.0000 ? ( 0.0242400*(exp(( - 0.0105200*STATES[0]))))/(1.00000+(exp(( - 0.137800*(STATES[0]+40.1400))))) : ( 0.600000*(exp(( 0.0570000*STATES[0]))))/(1.00000+(exp(( - 0.100000*(STATES[0]+32.0000))))));
ALGEBRAIC[39] = 1.00000/(ALGEBRAIC[18]+ALGEBRAIC[31]);
RATES[9] = (ALGEBRAIC[5] - STATES[9])/ALGEBRAIC[39];
ALGEBRAIC[6] = 1.00000/(1.00000+(exp(((- 8.00000 - STATES[0])/7.50000))));
ALGEBRAIC[19] = 1.40000/(1.00000+(exp(((- 35.0000 - STATES[0])/13.0000))))+0.250000;
ALGEBRAIC[32] = 1.40000/(1.00000+(exp(((STATES[0]+5.00000)/5.00000))));
ALGEBRAIC[40] = 1.00000/(1.00000+(exp(((50.0000 - STATES[0])/20.0000))));
ALGEBRAIC[42] =  1.00000*ALGEBRAIC[19]*ALGEBRAIC[32]+ALGEBRAIC[40];
RATES[11] = (ALGEBRAIC[6] - STATES[11])/ALGEBRAIC[42];
ALGEBRAIC[55] = (( (( CONSTANTS[21]*CONSTANTS[10])/(CONSTANTS[10]+CONSTANTS[22]))*STATES[2])/(STATES[2]+CONSTANTS[23]))/(1.00000+ 0.124500*(exp((( - 0.100000*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))))+ 0.0353000*(exp((( - STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1])))));
ALGEBRAIC[25] =  (( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*(log((CONSTANTS[11]/STATES[2])));
ALGEBRAIC[50] =  CONSTANTS[16]*(pow(STATES[7], 3.00000))*STATES[8]*STATES[9]*(STATES[0] - ALGEBRAIC[25]);
ALGEBRAIC[51] =  CONSTANTS[17]*(STATES[0] - ALGEBRAIC[25]);
ALGEBRAIC[56] = ( CONSTANTS[24]*( (exp((( CONSTANTS[27]*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))))*(pow(STATES[2], 3.00000))*CONSTANTS[12] -  (exp((( (CONSTANTS[27] - 1.00000)*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))))*(pow(CONSTANTS[11], 3.00000))*STATES[3]*CONSTANTS[26]))/( ((pow(CONSTANTS[29], 3.00000))+(pow(CONSTANTS[11], 3.00000)))*(CONSTANTS[28]+CONSTANTS[12])*(1.00000+ CONSTANTS[25]*(exp((( (CONSTANTS[27] - 1.00000)*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))))));
RATES[2] =  (( - 1.00000*(ALGEBRAIC[50]+ALGEBRAIC[51]+ 3.00000*ALGEBRAIC[55]+ 3.00000*ALGEBRAIC[56]))/( 1.00000*CONSTANTS[4]*CONSTANTS[2]))*CONSTANTS[3];
ALGEBRAIC[33] =  (( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*(log((CONSTANTS[10]/STATES[1])));
ALGEBRAIC[44] = 0.100000/(1.00000+(exp(( 0.0600000*((STATES[0] - ALGEBRAIC[33]) - 200.000)))));
ALGEBRAIC[45] = ( 3.00000*(exp(( 0.000200000*((STATES[0] - ALGEBRAIC[33])+100.000))))+(exp(( 0.100000*((STATES[0] - ALGEBRAIC[33]) - 10.0000)))))/(1.00000+(exp(( - 0.500000*(STATES[0] - ALGEBRAIC[33])))));
ALGEBRAIC[46] = ALGEBRAIC[44]/(ALGEBRAIC[44]+ALGEBRAIC[45]);
ALGEBRAIC[47] =  CONSTANTS[13]*ALGEBRAIC[46]* pow((CONSTANTS[10]/5.40000), 1.0 / 2)*(STATES[0] - ALGEBRAIC[33]);
ALGEBRAIC[54] =  CONSTANTS[20]*STATES[16]*STATES[15]*(STATES[0] - ALGEBRAIC[33]);
ALGEBRAIC[48] =  CONSTANTS[14]* pow((CONSTANTS[10]/5.40000), 1.0 / 2)*STATES[4]*STATES[5]*(STATES[0] - ALGEBRAIC[33]);
ALGEBRAIC[41] =  (( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*(log(((CONSTANTS[10]+ CONSTANTS[9]*CONSTANTS[11])/(STATES[1]+ CONSTANTS[9]*STATES[2]))));
ALGEBRAIC[49] =  CONSTANTS[15]*(pow(STATES[6], 2.00000))*(STATES[0] - ALGEBRAIC[41]);
ALGEBRAIC[52] = ( (( CONSTANTS[18]*STATES[11]*STATES[12]*STATES[13]*STATES[14]*4.00000*(STATES[0] - 15.0000)*(pow(CONSTANTS[2], 2.00000)))/( CONSTANTS[0]*CONSTANTS[1]))*( 0.250000*STATES[10]*(exp((( 2.00000*(STATES[0] - 15.0000)*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1])))) - CONSTANTS[12]))/((exp((( 2.00000*(STATES[0] - 15.0000)*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1])))) - 1.00000);
ALGEBRAIC[43] =  (( 0.500000*CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*(log((CONSTANTS[12]/STATES[3])));
ALGEBRAIC[53] =  CONSTANTS[19]*(STATES[0] - ALGEBRAIC[43]);
ALGEBRAIC[58] = ( CONSTANTS[32]*(STATES[0] - ALGEBRAIC[33]))/(1.00000+(exp(((25.0000 - STATES[0])/5.98000))));
ALGEBRAIC[57] = ( CONSTANTS[30]*STATES[3])/(STATES[3]+CONSTANTS[31]);
//ALGEBRAIC[12] = (VOI -  (floor((VOI/CONSTANTS[6])))*CONSTANTS[6]>=CONSTANTS[5]&&VOI -  (floor((VOI/CONSTANTS[6])))*CONSTANTS[6]<=CONSTANTS[5]+CONSTANTS[7] ? - CONSTANTS[8] : 0.00000);
RATES[0] = - (ALGEBRAIC[47]+ALGEBRAIC[54]+ALGEBRAIC[48]+ALGEBRAIC[49]+ALGEBRAIC[52]+ALGEBRAIC[55]+ALGEBRAIC[50]+ALGEBRAIC[51]+ALGEBRAIC[56]+ALGEBRAIC[53]+ALGEBRAIC[58]+ALGEBRAIC[57]+ALGEBRAIC[12]);
RATES[1] =  (( - 1.00000*((ALGEBRAIC[47]+ALGEBRAIC[54]+ALGEBRAIC[48]+ALGEBRAIC[49]+ALGEBRAIC[58]+ALGEBRAIC[12]) -  2.00000*ALGEBRAIC[55]))/( 1.00000*CONSTANTS[4]*CONSTANTS[2]))*CONSTANTS[3];
ALGEBRAIC[59] = CONSTANTS[44]/(1.00000+(pow(CONSTANTS[42], 2.00000))/(pow(STATES[3], 2.00000)));
ALGEBRAIC[60] =  CONSTANTS[43]*(STATES[17] - STATES[3]);
ALGEBRAIC[61] =  CONSTANTS[41]*(STATES[10] - STATES[3]);
ALGEBRAIC[63] = 1.00000/(1.00000+( CONSTANTS[45]*CONSTANTS[46])/(pow((STATES[3]+CONSTANTS[46]), 2.00000)));
RATES[3] =  ALGEBRAIC[63]*((( (ALGEBRAIC[60] - ALGEBRAIC[59])*CONSTANTS[51])/CONSTANTS[4]+ALGEBRAIC[61]) - ( 1.00000*((ALGEBRAIC[53]+ALGEBRAIC[57]) -  2.00000*ALGEBRAIC[56])*CONSTANTS[3])/( 2.00000*1.00000*CONSTANTS[4]*CONSTANTS[2]));
ALGEBRAIC[62] = CONSTANTS[38] - (CONSTANTS[38] - CONSTANTS[39])/(1.00000+(pow((CONSTANTS[37]/STATES[17]), 2.00000)));
ALGEBRAIC[65] =  CONSTANTS[34]*ALGEBRAIC[62];
RATES[18] =  - ALGEBRAIC[65]*STATES[10]*STATES[18]+ CONSTANTS[36]*(1.00000 - STATES[18]);
ALGEBRAIC[64] = CONSTANTS[33]/ALGEBRAIC[62];
ALGEBRAIC[66] = ( ALGEBRAIC[64]*(pow(STATES[10], 2.00000))*STATES[18])/(CONSTANTS[35]+ ALGEBRAIC[64]*(pow(STATES[10], 2.00000)));
ALGEBRAIC[67] =  CONSTANTS[40]*ALGEBRAIC[66]*(STATES[17] - STATES[10]);
ALGEBRAIC[68] = 1.00000/(1.00000+( CONSTANTS[47]*CONSTANTS[48])/(pow((STATES[17]+CONSTANTS[48]), 2.00000)));
RATES[17] =  ALGEBRAIC[68]*(ALGEBRAIC[59] - (ALGEBRAIC[67]+ALGEBRAIC[60]));
ALGEBRAIC[69] = 1.00000/(1.00000+( CONSTANTS[49]*CONSTANTS[50])/(pow((STATES[10]+CONSTANTS[50]), 2.00000)));
RATES[10] =  ALGEBRAIC[69]*((( - 1.00000*ALGEBRAIC[52]*CONSTANTS[3])/( 2.00000*1.00000*CONSTANTS[52]*CONSTANTS[2])+( ALGEBRAIC[67]*CONSTANTS[51])/CONSTANTS[52]) - ( ALGEBRAIC[61]*CONSTANTS[4])/CONSTANTS[52]);
}
}
