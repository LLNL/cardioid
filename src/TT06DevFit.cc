#include "TT06DevFit.hh"
#include <cmath>
#include <gsl/gsl_multifit.h>
#define C(i) (gsl_vector_get(c,(i)))
#define sigm(x)   ((x)/(1.0+(x)) )
static double cA[3],cB[3]; 
static double c1,c2,c3,c4,c5,c6,c7,c8,c9;
static double c10,c11,c12,c13,c14,c15,c16,c17,c18,c19;
static double c20,c21,c22,c23,c24,c25,c26,c27,c28,c29;
static double c30,c31,c32,c33,c34,c36,c40,c43,c44;
static double f1,f2,f3,f4,f5,f6,f7,f8,f9,f10; 
static  int fitL[30],fitM[30]; 
static double fitCoef[30][64]; 
void update_fCassGate(double dt, double Ca_SS, double *gate);
void  funcValue(int cellType, double Vm , double dt, double *gv); 
void  funcValueA(int cellType, double Vm , double dt, double *gv); 
void initState(double *STATES,int cellType)
{
if (cellType == 0) 
{
//STATES[Vmembrane] = -86.709;
STATES[K_i] = 138.4;
STATES[Na_i] = 10.355;
STATES[Ca_i] = 0.00013;
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
STATES[fCass_gate] = 0.9953;
STATES[Ca_ss] = 0.00036;
STATES[Ca_SR] = 3.715;
STATES[R_prime] = 0.9068;
}

if (cellType == 1) 
{
//STATES[Vmembrane] = -85.423;
STATES[K_i] = 138.52;
STATES[Na_i] = 10.132;
STATES[Ca_i] = 0.000153;
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
STATES[fCass_gate] = 0.9942;
STATES[Ca_ss] = 0.00042;
STATES[Ca_SR] = 4.272;
STATES[R_prime] = 0.8978;
}

if (cellType==2) 
{
//STATES[Vmembrane] = -85.23;
STATES[K_i] = 136.89;
STATES[Na_i] = 8.604;
STATES[Ca_i] = 0.000126;
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
STATES[fCass_gate] = 0.9953;
STATES[Ca_ss] = 0.00036;
STATES[Ca_SR] = 3.64;
STATES[R_prime] = 0.9073;
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
double get_c9() { return c9; }
double computeUpdates(double dt, double Vm, double* STATES, int cellType)
{
   double gv[30]; 
   funcValueA(cellType,Vm,dt,gv);
   double *fv= gv+24; 

#define logSeries(x) ((x)*(1.0-0.5*(x)))
   double dV0 = 1.0*Vm -c3*log(STATES[K_i]) -c5;
   double dV1 = 1.0*Vm -c3*log(STATES[Na_i]) -c4;
   //double dV2 = 1.0*Vm -c3*log(STATES[K_i]+c2*STATES[Na_i]) -c6;
   double dV3 = 2.0*Vm- c3*log(STATES[Ca_i]) - c8;
   //double ss = c2*(STATES[Na_i]/STATES[K_i]);
   //double dV2 =Vm-c3*log(STATES[K_i]) -c3*log(1+ss) -c6;
   //double dV2 =Vm-c3*log(STATES[K_i])-c3*logSeries(c2*STATES[Na_i]/STATES[K_i]) -c6;
   double dV2 =dV0-c3*logSeries(c2*STATES[Na_i]/STATES[K_i])+c5-c6;
   
   double xx  =  (3.0*exp(0.0002*dV0 + 0.02)+exp(0.1*dV0 - 1.0))/(1.0+exp( -0.5*dV0))*(10.0+10*exp(0.06*dV0 -12.0));
   double fdV0 =  c9*c12/(1.0+xx);
   
   double fdV3 =  c7*dV3;
   
   double x[8]; 
   x[0] = STATES[Na_i]*c25; 
   x[1] = STATES[Ca_i]*c26; 
   x[2] = SQ(STATES[Ca_i]*c27); 
   x[3] = SQ(STATES[Ca_i]*c28+c29); 
   x[4] = SQ(STATES[Ca_SR]*c30); 
   x[6] = SQ(STATES[Ca_SR]*c31+c32); 
   x[7] = SQ(STATES[Ca_ss]*c33+c34); 
   double sigm0 = sigm(x[0]); 
   double sigm1 = sigm(x[1]); 
   double sigm2 = sigm(x[2]); 
   double sigm3 = sigm(x[3]); 
   double sigm4 = sigm(x[4]); 
   double sigm6 = sigm(x[6]); 
   double sigm7 = sigm(x[7]); 
   
   
   double tmp5 =    fdV0 +      cB[cellType]*STATES[r_gate]*STATES[s_gate]+ c11*STATES[Xr1_gate]*STATES[Xr2_gate] +  fv[5];
   double tmp6 =  c20*CUBE(STATES[m_gate])*STATES[h_gate]*STATES[j_gate]+c21;
   double tmp7 =  cA[cellType]*SQ(STATES[Xs_gate]);
   double tmp8  = c18+c19*sigm4; //Sigm4
   
   double sigm5 =   SQ(STATES[Ca_ss])/(tmp8*c17 + SQ(STATES[Ca_ss]));
   double  tmp9 =    tmp8*STATES[Ca_ss]+c36; 
   
   double itmpA =  sigm0 * fv[0];                          //Sigm0
   double itmp0 = CUBE(STATES[Na_i])*fv[1]-STATES[Ca_i]*fv[2]; 
   double itmp1 =  itmp0-1.5*itmpA+tmp6*dV1; 
   double itmp2 = itmpA + tmp5*dV0+tmp7*dV2; 
   double itmp3 = STATES[d_gate]*STATES[f_gate]*STATES[f2_gate]*STATES[fCass_gate]*(fv[4]-STATES[Ca_ss]*fv[3]);
   double itmp4 = fdV3+c24*sigm1; 
   double itmp5=  c43*(STATES[Ca_i] - STATES[Ca_SR])+c44*sigm2;      
   double itmp6 = c23*(STATES[Ca_ss] - STATES[Ca_i]);
   double itmp7 = sigm5*STATES[R_prime]*(STATES[Ca_SR] - STATES[Ca_ss]);
   
   double dVdt  = itmp3+itmp1*c22+itmp2*c22-itmp4*c10;
   STATES[K_i]  += dt*(itmp2);
   STATES[Na_i]  += dt*(itmp1+2.0*itmp0);
   STATES[Ca_i]  += dt*(sigm3*(itmp4-itmp0+itmp6*c15-itmp5*c16));
   double Ca_SS = STATES[Ca_ss]; 
   STATES[Ca_ss] = Ca_SS + dt*(sigm7*(itmp6+itmp7*c14+itmp3*c13));    
   STATES[Ca_SR] += dt*(sigm6*(itmp5-c40*itmp7));
   STATES[R_prime] += dt*(c36 - tmp9*STATES[R_prime]);
   
   update_fCassGate(dt, Ca_SS, STATES+fCass_gate);

   for (int i=0;i<10;i++) STATES[i+Xr1_gate] = gv[2*i] + gv[2*i+1]*STATES[i+Xr1_gate]   ;
   int i  = 10 + (cellType != 0) ; 
   STATES[s_gate] = gv[2*i] + gv[2*i+1]*STATES[s_gate]   ;
   
   /*
   static FILE *file = NULL; 
   static double time = 0.0; 
   static int loop =0; 
   static double xmin[7],xmax[7]; 
   if (file == NULL) 
   {	file = fopen("info.data","w"); 
   //	for (int i=0;i<7;i++)  xmin[i]=xmax[i] = x[i]; 
   }
   if (loop % 10 ==0) fprintf(file,"%14.10f %e %24.14e %24.14e %24.14e %24.14e %24.14e %24.14e\n",time,RATES[0],RATES[K_i],RATES[Na_i],RATES[Ca_i],RATES[Ca_ss],RATES[Ca_SR],STATES[R_prime]); 
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
static void update_Xr1Gate(double dt, double Vm, double *gate,double *a, double *b)
{
   double mhu = 1.0/(1.0+(exp(((- 26.0 - Vm)/7.0))));
   double t1 = 450.0/(1.0+(exp(((- 45.0 - Vm)/10.0))));
   double t2 = 6.0/(1.0+(exp(((Vm+30.0)/11.5000))));
   double tau =  t1*t2;
   double rate = (mhu - *gate)/tau;
   *a =  mhu*dt/tau; 
   *b =  1.0-dt/tau; 
   *gate += rate*dt;
}
// Update Xr2 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
static void update_Xr2Gate(double dt, double Vm, double *gate, double *a, double *b)
{
   double mhu = 1.0/(1.0+(exp(((Vm+88.0)/24.0))));
   double t1 = 3.0/(1.0+(exp(((- 60.0 - Vm)/20.0))));
   double t2 = 1.12000/(1.0+(exp(((Vm - 60.0)/20.0))));
   double tau =  t1*t2;
   double rate = (mhu - *gate)/tau;
   *a =  mhu*dt/tau; 
   *b =  1.0-dt/tau; 
   *gate += rate*dt;
}
static void update_XsGate(double dt, double Vm, double *gate, double *a, double *b)
{
   double mhu = 1.0/(1.0+(exp(((- 5.0 - Vm)/14.0))));
   double t1 = 1400.00/ sqrt(1.0+exp((5.0 - Vm)/6.0));
   double t2 = 1.0/(1.0+(exp(((Vm - 35.0)/15.0))));
   double tau  =  t1*t2+80.0;
   double rate = (mhu - *gate)/tau;
   *a =  mhu*dt/tau; 
   *b =  1.0-dt/tau; 
   *gate += rate*dt;
}
static void update_mGate(double dt, double Vm, double *gate,double *a, double *b)
{
   double mhu = 1.0/SQ(1.0+exp((-56.8600 - Vm)/9.03000));
   double t1  = 1.0/(1.0+(exp(((- 60.0 - Vm)/5.0))));
   double t2  = 0.100000/(1.0+(exp(((Vm+35.0)/5.0))))+0.100000/(1.0+(exp(((Vm - 50.0)/200.0))));
   double tau =  t1*t2;
   //double rate = (mhu - *gate)/tau;
   *a = mhu*(1-exp(-dt/tau)); 
   *b = exp(-dt/tau); 
   *gate = mhu - (mhu -*gate) *exp(-dt/tau);
}
static void update_hGate(double dt, double Vm, double *gate,double *a, double *b)
{
   double mhu = 1.0/SQ((1.0+(exp(((Vm+71.5500)/7.43000)))));
   double t1  = (Vm<- 40.0 ?  0.0570000*(exp((- (Vm+80.0)/6.80000))) : 0.0);
   double t2  = (Vm<- 40.0 ?  2.70000*(exp(( 0.0790000*Vm)))+ 310000.*(exp(( 0.348500*Vm))) : 0.770000/( 0.130000*(1.0+(exp(((Vm+10.6600)/- 11.1000))))));
   double tau = 1.0/(t1+t2);
   *a =  mhu*dt/tau; 
   *b =  1.0-dt/tau; 
   double rate = (mhu - *gate)/tau;
   *gate += rate*dt;
}
static void update_jGate(double dt, double Vm, double *gate,double *a, double *b)
{
   double mhu = 1.0/SQ((1.0+(exp(((Vm+71.5500)/7.43000)))));
   double t1  = (Vm < -40.0 ? (( ( - 25428.0*(exp(( 0.244400*Vm))) -  6.94800e-06*(exp(( - 0.0439100*Vm))))*(Vm+37.7800))/1.0)/(1.0+(exp(( 0.311000*(Vm+79.2300))))) : 0.0);
   double t2 = (Vm < -40.0 ? ( 0.0242400*(exp(( - 0.0105200*Vm))))/(1.0+(exp(( - 0.137800*(Vm+40.1400))))) : ( 0.600000*(exp(( 0.0570000*Vm))))/(1.0+(exp(( - 0.100000*(Vm+32.0))))));
   double tau  = 1.0/(t1+t2);
   *a =  mhu*dt/tau; 
   *b =  1.0-dt/tau; 
   double rate = (mhu - *gate)/tau;
   *gate += rate*dt;
}
static void update_rGate(double dt, double Vm, double *gate,double *a, double *b)
{
   double mhu = 1.0/(1.0+(exp(((20.0 - Vm)/6.0))));
   double tau =  9.50000*(exp((- SQ((Vm+40.0)))/1800.00))+0.800000;
   *a =  mhu*dt/tau; 
   *b =  1.0-dt/tau; 
   double rate = (mhu - *gate)/tau;
   *gate += rate*dt;
}
static void update_dGate(double dt, double Vm, double *gate,double *a, double *b)
{
   double mhu = 1.0/(1.0+(exp(((- 8.0 - Vm)/7.50000))));
   double t1  = 1.40000/(1.0+(exp(((- 35.0 - Vm)/13.0))))+0.250000;
   double t2 = 1.40000/(1.0+(exp(((Vm+5.0)/5.0))));
   double t3 = 1.0/(1.0+(exp(((50.0 - Vm)/20.0))));
   double tau =  t1*t2+t3;
   *a =  mhu*dt/tau; 
   *b =  1.0-dt/tau; 
   double rate = (mhu - *gate)/tau;
   *gate += rate*dt;
}
static void update_fGate(double dt, double Vm, double *gate,double *a, double *b)
{
   double mhu = 1.0/(1.0+(exp(((Vm+20.0)/7.0))));
   double tau =  1102.50*(exp((- SQ(Vm+27.0)/225.0)))+200.0/(1.0+(exp(((13.0 - Vm)/10.0))))+180.0/(1.0+(exp(((Vm+30.0)/10.0))))+20.0;
   *a =  mhu*dt/tau; 
   *b =  1.0-dt/tau; 
   double rate = (mhu - *gate)/tau;
   *gate += rate*dt;
}
static void update_f2Gate(double dt, double Vm, double *gate,double *a, double *b)
{
   double mhu = 0.670000/(1.0+(exp(((Vm+35.0)/7.0))))+0.330000;
   double tau =  562.0*exp(-SQ((Vm+27.0))/240.0)+31.0/(1.0+(exp(((25.0 - Vm)/10.0))))+80.0/(1.0+(exp(((Vm+30.0)/10.0))));
   *a =  mhu*dt/tau; 
   *b =  1.0-dt/tau; 
   double rate = (mhu - *gate)/tau;
   *gate += rate*dt;
}
static void update_sGate(double dt, double Vm, double *gate, int cellType,double *a, double *b)
{

   if (cellType==0) 
   {
      double mhu = 1.0/(1.0+(exp(((Vm+28.0)/5.0))));
      double tau =  1000.0*(exp((- (pow((Vm+67.0), 2.0))/1000.0)))+8.0;
   *a =  mhu*dt/tau; 
   *b =  1.0-dt/tau; 
      double rate = (mhu - *gate)/tau;
      *gate += rate*dt;
   }
   else
   {
      double mhu = 1.0/(1.0+(exp(((Vm+20.0)/5.0))));
      double tau =  85.0*(exp((- (pow((Vm+45.0), 2.0))/320.0)))+5.0/(1.0+(exp(((Vm - 20.0)/5.0))))+3.0;
   *a =  mhu*dt/tau; 
   *b =  1.0-dt/tau; 
      double rate = (mhu - *gate)/tau;
      *gate += rate*dt;
   }
}
void update_fCassGate(double dt, double Ca_SS, double *gate)
{
   double t1 = 1.0/(1.0+SQ(20*Ca_SS)); 
   double mhu = 0.600000*t1+0.4000000;
   double tau =    80.0*t1+2.0;
   double rate = (mhu - *gate)/tau;
   *gate += rate*dt;
}
void extra(double Vm, double *fv)
{
double expV1 = exp(-f1*Vm); 
double expV2 = exp(f1*(30-2.0*Vm)); 
double expV5 = exp(f5*Vm); 
double expV = expV5*expV1; 
double fva = 1.0/(f3+ f4*expV);
fv[0]= f2/(1.0+0.1245*exp(-0.1*Vm*f1)+0.0353*expV1);
fv[1] = expV5*fva ; 
fv[2] = expV* fva*f6 ; 
fv[3] = (f7*Vm-f8)/(1.0-expV2);
fv[4] = fv[3] * expV2 *f9;
fv[5] = f10/(1.0+exp((25.0 - Vm)/5.98));
}

int costFunc(int l, int m)
{
   int cost = (m-1)+(l-1)+4;
   if ( l==1 ) cost = (m-1) ; 
	return cost; 
}
double fit(int l,int m,double  *a,double x, double *s1, double *s2) 
{
   double sum1=0; 
   double sum2=0;    
   int k = m+l-1; 
   for (int j=m-1;j>=0;j--)sum1 =  a[j] + x*sum1; 
   if (l<2) return sum1; 
   for (int j=k;j>=m;j--)  sum2 =  a[j] + x*sum2; 
   *s1=sum1; 
   *s2=sum2; 
   return sum1/sum2; 
}
void  funcValue(int cellType, double Vm , double dt, double *gv)
{
   double x0; 
   int i=0;
   x0 = 0; update_Xr1Gate(dt, Vm, &x0, gv+i, gv+i+1);    i+=2; 
   x0 = 0; update_Xr2Gate(dt, Vm, &x0, gv+i, gv+i+1);    i+=2; 
   x0 = 0; update_XsGate(dt, Vm, &x0, gv+i, gv+i+1);    i+=2; 
   x0 = 0; update_mGate(dt, Vm, &x0, gv+i, gv+i+1);    i+=2; 
   x0 = 0; update_hGate(dt, Vm, &x0, gv+i, gv+i+1);    i+=2; 
   x0 = 0; update_jGate(dt, Vm, &x0, gv+i, gv+i+1);    i+=2; 
   x0 = 0; update_rGate(dt, Vm, &x0, gv+i, gv+i+1);    i+=2; 
   x0 = 0; update_dGate(dt, Vm, &x0, gv+i, gv+i+1);    i+=2; 
   x0 = 0; update_fGate(dt, Vm, &x0, gv+i, gv+i+1);    i+=2; 
   x0 = 0; update_f2Gate(dt, Vm, &x0, gv+i, gv+i+1);    i+=2; 
   x0 = 0; update_sGate(dt, Vm, &x0, 0 ,gv+i, gv+i+1);    i+=2; 
   x0 = 0; update_sGate(dt, Vm, &x0, 1,gv+i, gv+i+1);    i+=2; 
   extra(Vm,gv+i); 
}
void  funcValueA(int cellType, double Vm , double dt, double *gv)
{
   double s1,s2; 
   for (int i=0;i<30;i++) gv[i] = fit(fitL[i],fitM[i],fitCoef[i],Vm,&s1,&s2); 
}


   

void calcError(int funcNumber, int l,int m,double *a,int n,double *x,double *y,double *errMax, double *errRMS,double *error)
{
   double eMax = 0.0; 
   double err2 = 0.0; 
   double s1,s2; 
   for (int i = 0; i<n;i++) 
   {
     double f = fit(l,m,a,x[i],&s1,&s2); 
     double err = fabs(f-y[i]); 
     error[i]=err; 
     if (err > eMax) eMax = err; 
     err2 += err*err; 
   }
   *errMax = eMax; 
   *errRMS = sqrt(err2/n); 
}
double  calcQ(int funcNumber, double x)
{
   double dt = 0.02; 
   double gv[30]; 
   funcValue(0,x,dt,gv);
   return gv[funcNumber]; 
}
void  pade(int l, int m, int n, double *x, double *y, double *a )
{
   gsl_matrix *X=NULL, *cov=NULL;
   gsl_vector *yy=NULL, *w=NULL, *c=NULL;
   int k = m+l-1 ; 
      
    X = gsl_matrix_alloc(n, k);
    yy= gsl_vector_alloc(n);
    w = gsl_vector_alloc(n);
        
     c = gsl_vector_alloc(k);
     cov = gsl_matrix_alloc(k, k);
         
     for (int i = 0; i < n; i++)
     {
           double xj=1.0; 
           for (int j=0;j<m;j++) { gsl_matrix_set (X, i, j, xj); xj *= x[i]; }
           xj=y[i]*x[i]; 
           for (int j=m;j<k;j++) { gsl_matrix_set (X, i, j, xj); xj *= x[i]; }
           gsl_vector_set (yy, i, y[i]);
           gsl_vector_set (w, i, 1.0);
      }
      for (int i = 0; i < n; i++) gsl_vector_set(w, i, 1.0);
      double chisq; 
      {
          gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, k);
          gsl_multifit_wlinear (X, w, yy, c, cov, &chisq, work);
          gsl_multifit_linear_free (work);
       }
       for (int j=0;j<m;j++) a[j] = C(j);
       a[m]=1.0; 
       for (int j=m;j<k;j++) a[j+1] =-C(j);
      gsl_matrix_free (X);
      gsl_vector_free (yy);
      gsl_vector_free (w);
      gsl_vector_free (c);
      gsl_matrix_free (cov);
}
void makeFunctionTable(int funcNumber, double x0, double x1, int n, double *x, double *y, double *yMin, double *yMax, double *dydxMax) 
{
   double deltaX = (x1-x0)/n; 
   n++; 
   double ymin=0.0,ymax=0.0,dydxmax=0.0; 
   for (int i = 0; i < n; i++)
   {
      x[i] = x0 + deltaX * i ; 
      y[i]=calcQ(funcNumber,x[i]); 
      double dydx = (y[i]-y[i-1])/(x[i]-x[i-1]); 
      if (dydx > dydxmax || i==0) dydxmax = dydx;   
      if ( y[i] < ymin   || i==0) ymin    =  y[i]; 
      if ( y[i] > ymax   || i==0) ymax    =  y[i]; 
   }
   *yMin = ymin; 
   *yMax = ymax; 
   *dydxMax = dydxmax; 
}
void errorInfo(int funcNumber, double tol, int l, int m, int n, double *a, double *x, double *y, double yMin, double yMax, double dydxMax) 
{
   double dy = yMax-yMin; 
   double tableRes = 0.1; // millivolts
   double tableErrMax = 0.5*tableRes*dydxMax; 
   double errMax,errRMS,err[n+1]; 
   calcError(funcNumber,l,m,a,n,x,y,&errMax,&errRMS,err);
   fprintf(stderr,"# %4d: %6d %6d %6d ",funcNumber,costFunc(l,m),l,m); fflush(stdout); 
   fprintf(stderr,"%10.2e %10.2e %10.2e ",tol,errMax,errMax/dy); 
   fprintf(stderr,"%10.2e %10.2e %10.2e\n",errRMS,errRMS/dy,tableErrMax/dy); 
   char filename[256]; 
   sprintf(filename,"func_tt06.%d",funcNumber); 
   FILE *file = fopen(filename,"w"); 
   for (int i=0;i<n;i++) 
   {
      double s1,s2; 
      double f = fit(l,m,a,x[i],&s1,&s2); 
      fprintf(file,"%e %e %e %e\n",x[i],y[i],f,err[i]); fflush(file); 
   }
}
     
void  fitVoltageFunc(int funcNumber, double errTol, double *x, double *y, int n, int lMax, int mMax, int *lMin, int *mMin, double *amin)
{
   double err[n]; 
   int lmin=0.0,mmin=0.0; 
   for (int kk=4;kk<128;kk++) 
   {
   	double errMin = -1.0; 
        for (int l=1;l<=lMax;l+=1)
        for (int m=1;m<=mMax;m+=1)
        {
           //if (l != m  ) continue; 
	   //if (l > 1 ) continue ; 
           if (costFunc(l,m) != kk) continue; 
           double errMax, errRMS; 
           double a[l+m]; 
	   pade(l,m,n,x, y,a);
           calcError(funcNumber,l,m,a,n,x,y,&errMax,&errRMS,err);
           if (errMax < errMin || errMin==-1.0) 
           {
              lmin=l;
              mmin=m;
              for (int j=0;j<l+m;j++) amin[j] = a[j];
              errMin = errMax; 
           }
       }
       if (errMin  < errTol && errMin >= 0.0 ) break ; 
   }
   *lMin = lmin; 
   *mMin = mmin; 
}
void  Approx (int n, int maxOrder, double tol)
{
   
   int lMax = maxOrder; 
   int mMax = maxOrder; 
   

   double x0 = -90; 
   double x1 =  40; 
   double yMin,yMax, dydxMax; 
   double x[n+1],y[n+1]; 
   int cost =0; 
   for (int p =0;p<30;p++) 
   { 
   makeFunctionTable(p, x0,x1,n, x, y, &yMin, &yMax, &dydxMax) ;
   fitVoltageFunc(p,tol*(yMax-yMin),x,y,n,lMax,mMax,&fitL[p],&fitM[p],fitCoef[p]);
   errorInfo(p,tol, fitL[p], fitM[p], n, fitCoef[p], x, y, yMin, yMax, dydxMax) ;
   cost += costFunc(fitL[p],fitM[p]); 
   }
   fprintf(stderr,"Total Cost = %d\n", cost); 
}
