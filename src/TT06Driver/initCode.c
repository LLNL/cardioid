#include <math.h>
#include "TT06Func.h" 
#define NFIT 33 
static double P_NaK,g_Ks[3],g_to[3]; 
static double g_NaL =0.0; 
//static double cA[3],cB[3]; 
static double c1,c2,c3,c4,c5,c6,c7,c8,c9;
static double c10,c11,c12,c13,c14,c15,c16,c17,c18,c19;
static double c20,c21,c22,c23,c24,c25,c26,c27,c28,c29;
static double c30,c31,c32,c33,c34,c36,c40,c43,c44;
static double f1,f2,f3,f4,f5,f6,f7,f7a,f9,f9a,f10; 
void (*fv05Func)(double Vm, double *fv);
double (*fv6Func)(double dv);

// JAG static double c[] ={0.5*M_LN2,2.0,2/3.0,2/5.0,2.0/7.0,2.0/9.0,2.0/11.0,2/13.0}; // c[0] = 0.5*log(2)


static double p[32];
static double d[32];

/*
static double Logc0 = 0.5*M_LN2;
static double Logc1 = 2.0;
static double Logc2 = 2.0/3.0;
static double Logc3 = 2.0/5.0;
static double Logc4 = 2.0/7.0;
static double Logc5 = 2.0/9.0;
static double Logc6 = 2.0/11.0;
static double Logc7 = 2.0/13.0;

static   double v_array[40]__attribute__((aligned(32)));
static   double arrayM_SQRT1_2[4]__attribute__((aligned(32)));  
static   double arrayM_LN2[4]__attribute__((aligned(32)));  
static   double array1[4]__attribute__((aligned(32)));  

void initArray()
{
  for (int ii=0; ii<4; ++ii) {
    v_array[ii]    = Logc0;
    v_array[ii+4]  = Logc1;
    v_array[ii+8]  = Logc2;
    v_array[ii+12] = Logc3;
    v_array[ii+16] = Logc4;
    v_array[ii+20] = Logc5;
    v_array[ii+24] = Logc6;
    v_array[ii+28] = Logc7;
    arrayM_SQRT1_2[ii] = M_SQRT1_2;
    arrayM_LN2[ii]     = M_LN2;
    array1[ii]         = 1.0;
  }
}

*/

void fastLogInitJAG()
{
    for (int n=0;n<32;n++) 
    {
       p[n] =  pow(0.5,n-16); 
       d[n] =  (n-16)*M_LN2; 
    }
}






#include <math.h>
#include <stdio.h>
#include <sys/param.h>
#ifdef BYTE_ORDER
# if BYTE_ORDER ==LITTLE_ENDIAN
#define LeadingBits 1
# else
#  if BYTE_ORDER == BIG_ENDIAN
#define LeadingBits 0
#  else
     Error: unknown byte order!
#  endif
# endif
#endif 
//#define fastLog(x) x
inline double fastLog(double x)
{
double c[] ={0.5*M_LN2,2.0,2/3.0,2/5.0,2.0/7.0,2.0/9.0,2.0/11.0,2/13.0}; // c[0] = 0.5*log(2); 
       union {unsigned int u[2]; double d ;} t;
       t.d=x; 
       int n = (t.u[LeadingBits]>>20)-1023; 
       t.u[LeadingBits] -= n<<20;
       double m =  t.d;  
       double z = (M_SQRT1_2*m-1.0)/(M_SQRT1_2*m+1.0);
       double z2 = z*z; 
       double f = n*M_LN2 +c[0] + z*(c[1] + z2*(c[2] + z2*(c[3]+z2*(c[4]+z2*(c[5]+z2*(c[6]+z2*c[7]))))));
       return f ; 
}



char *getStateName(int index)
{

char *stateNames [] = { "Ca_i"  , "K_i"    , "Na_i"    , "Ca_ss"   , "Ca_SR"  , "R_prime", "fCass" , "m_gate",          
  "h_gate", "j_gate" , "Xr1_gate", "Xr2_gate", "Xs_gate", "r_gate" , "d_gate", "f_gate", "f2_gate", "jL_gate" , "s_gate"  } ; 
return stateNames[index]; 
}





double initState(double *states,double *gates, int cellType)
{
double Vm; 
if (cellType == 0) 
{
Vm = -86.709;
states[dVK_i] = 138.4/c9+Vm;
//states[K_i] = 138.4;
states[Na_i] = 10.355;
states[Ca_i] = 0.00013;
states[Ca_ss] = 0.00036;
states[Ca_SR] = 3.715;
states[R_prime] = 0.9068;
states[fCass] = 0.9953;
states[Xr1_gate] = 0.00448;
states[Xr2_gate] = 0.476;
states[Xs_gate] = 0.0087;
states[m_gate] = 0.00155;
states[h_gate] = 0.7573;
states[j_gate] = 0.7225;
states[r_gate] = 2.235e-8;
states[d_gate] = 3.164e-5;
states[f_gate] = 0.8009;
states[f2_gate]= 0.9778;
states[s_gate] = 0.3212;
states[jL_gate] = 0.066;
}

if (cellType == 1) 
{
Vm = -85.423;
states[dVK_i] = 138.52/c9+Vm;
//states[K_i] = 138.52;
states[Na_i] = 10.132;
states[Ca_i] = 0.000153;
states[Ca_ss] = 0.00042;
states[Ca_SR] = 4.272;
states[R_prime] = 0.8978;
states[fCass] = 0.9942;
states[Xr1_gate] = 0.0165;
states[Xr2_gate] = 0.473;
states[Xs_gate] = 0.0174;
states[m_gate] = 0.00165;
states[h_gate] = 0.749;
states[j_gate] = 0.6788;
states[r_gate] = 2.347e-8;
states[d_gate] = 3.288e-5;
states[f_gate] = 0.7026;
states[f2_gate] = 0.9526;
states[s_gate] = 0.999998;
states[jL_gate] = 0.066;
}

if (cellType==2) 
{
Vm  = -85.23;
states[dVK_i] = 136.89/c9+Vm;
//states[K_i] = 136.89;
states[Na_i] = 8.604;
states[Ca_i] = 0.000126;
states[Ca_ss] = 0.00036;
states[Ca_SR] = 3.64;
states[R_prime] = 0.9073;
states[fCass] = 0.9953;
states[Xr1_gate] = 0.00621;
states[Xr2_gate] = 0.4712;
states[Xs_gate] = 0.0095;
states[m_gate] = 0.00172;
states[h_gate] = 0.7444;
states[j_gate] = 0.7045;
states[r_gate] = 2.42e-8;
states[d_gate] = 3.373e-5;
states[f_gate] = 0.7888;
states[f2_gate] = 0.9755;
states[s_gate] = 0.999998;
states[jL_gate] = 0.066;
}
if (cellType == 3) 
{
Vm = -85.423;
states[dVK_i] = 138.52/c9+Vm;
//states[K_i] = 138.52;
states[Na_i] = 10.132;
states[Ca_i] = 0.000153;
states[Ca_ss] = 0.00042;
states[Ca_SR] = 4.272;
states[R_prime] = 0.8978;
states[fCass] = 0.9942;
states[Xr1_gate] = 0.0165;
states[Xr2_gate] = 0.473;
states[Xs_gate] = 0.0174;
states[m_gate] = 0.00165;
states[h_gate] = 0.749;
states[j_gate] = 0.6788;
states[r_gate] = 2.347e-8;
states[d_gate] = 3.288e-5;
states[f_gate] = 0.7026;
states[f2_gate] = 0.9526;
states[s_gate] = 0.999998;
states[jL_gate] = 0.066;
}
return Vm; 
}
double SP[40]__attribute__((aligned(32)));
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
   c7 =  (0.50*cnst[19]);
   c11= cnst[14]*sqrt(cnst[10]/5.4);
   c12= cnst[13]*sqrt(cnst[10]/5.4);
   c13= cnst[3]/(2.0*cnst[52]*cnst[2]*c9);
   c14 = cnst[51]*cnst[40]/(cnst[52]*c9); 
   c15 = -cnst[52]/cnst[4]; 
   c16 = cnst[51]/cnst[4];
   c17 = cnst[35]/(cnst[33]*cnst[34])*c9;
   c18 = cnst[34]*cnst[38]/c9 ;
   c19  = -cnst[34]*(cnst[38]-cnst[39])/c9; 
   c20  = cnst[16]; 
   c21  = cnst[17]; 
   c22  = 1/c9; 
   c23  = cnst[41]/(c15*c9); 
   //c24  = cnst[30]/c10; 
   c24  =  0.5*cnst[30]; 
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
   c36  =  cnst[36]/c9; 
   c40  =  cnst[40]/c9; 
   c43  =  cnst[43]/c9; 
   c44  =  cnst[44]/c9; 
   
   P_NaK= cnst[21]; 
   g_NaL =0.0; 
   
   f1 = c1; 
   f2 =  -2.0*cnst[10]/(cnst[10]+cnst[22]);
   f3 =  ((CUBE(cnst[29])+CUBE(cnst[11]))*(cnst[28]+cnst[12]))/(cnst[24]*cnst[12]); 
   f4 =  f3*cnst[25]; 
   f5 =  cnst[27]*f1; 
   f6 =  (CUBE(cnst[11])*cnst[26]/cnst[12]);
   f7 = cnst[18]*cnst[2]*f1;
   f7a = 0.5*cnst[18]*cnst[2];
   f9  = 4.0*cnst[12];
   f9a = 4.0*cnst[12]*f7a;
   f10 = cnst[32];
 SP[0] = c26 ;
 SP[1] = c27 ;
 SP[2] = c28 ;
 SP[3] = c29 ;
 SP[4] = c8 ;
 SP[5] = c7 ;
 SP[6] = c24 ;
 SP[7] = c43 ;
 SP[8] = c44 ;
 SP[9] = c23 ;
 SP[10] = c15 ;
 SP[11] = c16 ;
 SP[12] = c9 ;
 SP[13] = c25 ;
 SP[14] = c3 ;
 SP[15] = c5 ;
 SP[16] = c4 ;
 SP[17] = c2 ;
 SP[18] = c6 ;
 SP[19] = c11 ;
 SP[20] = c20 ;
 SP[21] = c21 ;
 SP[22] = c22 ;
 SP[23] = c30 ;
 SP[24] = c31 ;
 SP[25] = c32 ;
 SP[26] = c33 ;
 SP[27] = c34 ;
 SP[28] = c19 ;
 SP[29] = c18 ;
 SP[30] = c36 ;
 SP[31] = c17 ;
 SP[32] = 20.0 ;
 SP[33] = 0.6 ;
 SP[34] = 0.4 ;
 SP[35] = 80.0 ;
 SP[36] = c14 ;
 SP[37] = c13 ;
 SP[38] = c40 ;
 SP[39] = c36 ;
}

