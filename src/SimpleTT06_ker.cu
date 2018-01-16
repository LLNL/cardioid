
//#include "jitifyDecl.hpp"

namespace SimpleTT06 {
   
enum StateOffset {
   f2Gate_off,
   fGate_off,
   dGate_off,
   mGate_off,
   jGate_off,
   hGate_off,
   rGate_off,
   sGate_off,
   Xr1Gate_off,
   Xr2Gate_off,
   XsGate_off,
   jLGate_off,
   Na_i_off,
   Ca_i_off,
   Ca_ss_off,
   Ca_sr_off,
   fCass_off,
   dVK_i_off,
   R_prime_off,
   NUMSTATES
};

#define CUBE(x) ((x)*(x)*(x))
#define SQ(x) ((x)*(x))
#define sigm(x)   ((x)/(1+(x)))
#define logSeries(x)    (log(1+(x)) )

__global__ void actualCalcCuda(const double dt, const int nCells_, const double Vm[], const double iStim[], double dVm[], double* stateData)
{   
   int ii = blockIdx.x*blockDim.x + threadIdx.x;
   if (ii > nCells_) { return; }
   //printf("CPU state = %p\n", &state);
   //printf("CPU state->Ca_i = %p\n", state->Ca_i);
   //printf("CPU state->nCells = %d\n", state->nCells);
   //printf("CPU istim[0] = %g\n", iStim[0]);

      //printf("CPU state = %p\n", &state);
      //printf("CPU state->Ca_i = %p\n", state->Ca_i);
      //printf("CPU state->nCells = %d\n", state->nCells);
      //printf("GPU istim[0] = %g\n", iStim[0]);


      //set Vm
      const double thisVm = Vm[ii];
      const double istim = iStim[ii];

      //set all state variables
      //EDIT_STATE
      const double f2Gate=stateData[f2Gate_off*nCells_+ii];
      const double fGate=stateData[fGate_off*nCells_+ii];
      const double dGate=stateData[dGate_off*nCells_+ii];
      const double mGate=stateData[mGate_off*nCells_+ii];
      const double jGate=stateData[jGate_off*nCells_+ii];
      const double hGate=stateData[hGate_off*nCells_+ii];
      const double rGate=stateData[rGate_off*nCells_+ii];
      const double sGate=stateData[sGate_off*nCells_+ii];
      const double Xr1Gate=stateData[Xr1Gate_off*nCells_+ii];
      const double Xr2Gate=stateData[Xr2Gate_off*nCells_+ii];
      const double XsGate=stateData[XsGate_off*nCells_+ii];
      const double jLGate=stateData[jLGate_off*nCells_+ii];
      const double _Na_i=stateData[Na_i_off*nCells_+ii];
      const double _Ca_i=stateData[Ca_i_off*nCells_+ii];
      const double _Ca_ss=stateData[Ca_ss_off*nCells_+ii];
      const double _Ca_SR=stateData[Ca_sr_off*nCells_+ii];
      const double _fCass=stateData[fCass_off*nCells_+ii];
      const double _dVK_i=stateData[dVK_i_off*nCells_+ii];
      const double _R_prime=stateData[R_prime_off*nCells_+ii];

      //set per-cell flags
      //EDIT_PERCELL_FLAGS
      
      //set per-cell parameters
      //EDIT_PERCELL_PARAMETERS

      const double pcnst_0 = 8314.472;
      const double pcnst_1 = 310;
      const double pcnst_2 = 96485.3415;
      const double pcnst_9 = 0.03;
      const double pcnst_10 = 5.4;
      const double pcnst_11 = 140;
      const double pcnst_12 = 2;
      const double c1 = pcnst_2/(pcnst_0*pcnst_1);
      const double c2 = pcnst_9;
      const double c3 = -1/c1;

      //const double c4 = -c3*log(pcnst_11);
      //const double c5 = -c3*log(pcnst_10);
      //const double c6 = -c3*log(pcnst_10+pcnst_9*pcnst_11);
      //const double c8 = -0.5*c3*log(pcnst_12);
      //if (0) {
      //  printf("%.17g %.17g %.17g %.17g\n", c4,c5,c6,c8);
      //}

      const double c4 = 132.00985294338352;
      const double c5 = 45.050058022436289;
      const double c6 = 60.420198081560486;
      const double c8 = 9.2582839417106104;
      const double P_NaK=2.724;
      const double g_Ks = 0.392;
      const double g_Kr = 0.153;
      const double g_to = 0.073;
      const double g_NaL=0.0;

      const double pcnst_3 = 0.185;
      const double pcnst_4 = 0.016404;
      const double pcnst_5 = 10;
      const double pcnst_6 = 1000;
      const double pcnst_7 = 1;
      const double pcnst_8 = 52;
      const double pcnst_13 = 5.405;
      const double pcnst_14 = 0.153;
      const double pcnst_16 = 14.838;
      const double pcnst_17 = 0.00029;
      const double pcnst_18 = 0.0000398;
      const double pcnst_19 = 0.000592;
      const double pcnst_21 = 2.724;
      const double pcnst_22 = 1;
      const double pcnst_23 = 40;
      const double pcnst_24 = 1000;
      const double pcnst_25 = 0.1;
      const double pcnst_26 = 2.5;
      const double pcnst_27 = 0.35;
      const double pcnst_28 = 1.38;
      const double pcnst_29 = 87.5;
      const double pcnst_30 = 0.1238;
      const double pcnst_31 = 0.0005;
      const double pcnst_32 = 0.0146;
      const double pcnst_33 = 0.15;
      const double pcnst_34 = 0.045;
      const double pcnst_35 = 0.06;
      const double pcnst_36 = 0.005;
      const double pcnst_37 = 1.5;
      const double pcnst_38 = 2.5;
      const double pcnst_39 = 1;
      const double pcnst_40 = 0.102;
      const double pcnst_41 = 0.0038;
      const double pcnst_42 = 0.00025;
      const double pcnst_43 = 0.00036;
      const double pcnst_44 = 0.006375;
      const double pcnst_45 = 0.2;
      const double pcnst_46 = 0.001;
      const double pcnst_47 = 10;
      const double pcnst_48 = 0.3;
      const double pcnst_49 = 0.4;
      const double pcnst_50 = 0.00025;
      const double pcnst_51 = 0.001094;
      const double pcnst_52 = 0.00005468;

      const double c9 = -pcnst_3/(pcnst_4*pcnst_2);


      const double c7 =  pcnst_19;
      //c11= pcsnt[14]*sqrt(pcnst_10/5.4);
      //const double c11= sqrt(pcnst_10/5.4);
      const double c11 = 1;
      //const double c12= pcnst_13*sqrt(pcnst_10/5.4);
      const double c12 = 5.4050000000000002;
      const double c13= pcnst_3/(2.0*pcnst_52*pcnst_2*c9);
      const double c14 = pcnst_51/pcnst_52;
      const double c15 = -pcnst_52/pcnst_4;
      const double c16 = pcnst_51/pcnst_4;
      const double c17 = pcnst_35/(pcnst_33*pcnst_34)*c9;
      const double c18 = pcnst_34*pcnst_38/c9 ;
      const double c19  = -pcnst_34*(pcnst_38-pcnst_39)/c9;
      const double c20  = pcnst_16;
      const double c21  = pcnst_17;
      const double c22  = 1/c9;
      const double c23  = pcnst_41/(c15*c9);
      const double c24  =  pcnst_30;
      const double c25  =  1.0/pcnst_23;
      const double c26  =  1.0/pcnst_31;
      const double c27  =  1.0/pcnst_42;
      //const double c28  =  1.0/sqrt(pcnst_45*pcnst_46);
      const double c28 = 70.710678118654755;
      const double c29  =  pcnst_46*c28;
      const double c30  =  1.0/pcnst_37;
      //const double c31  =  1.0/sqrt(pcnst_47*pcnst_48);
      const double c31 = 0.57735026918962584;
      const double c32  =  pcnst_48*c31;
      //const double c33  =  1.0/sqrt(pcnst_49*pcnst_50);
      const double c33 = 100;
      const double c34  =  pcnst_50*c33;
      const double c36  =  pcnst_36/c9;
      const double c40  =  pcnst_40/c9;
      const double c43  =  pcnst_43/c9;
      const double c44  =  pcnst_44/c9;

      const double f1 = c1;
      const double f2 =  -2.0*pcnst_10/(pcnst_10+pcnst_22);
      const double f3 =  ((CUBE(pcnst_29)+CUBE(pcnst_11))*(pcnst_28+pcnst_12))/(pcnst_24*pcnst_12);
      const double f4 =  f3*pcnst_25;
   
      const double f5 =  pcnst_27*f1;
      const double f6 =  (CUBE(pcnst_11)*pcnst_26/pcnst_12);
      const double f7 = pcnst_18*pcnst_2*f1;
      const double f7a = 0.5*pcnst_18*pcnst_2;
      const double f9  = 4.0*pcnst_12;
      const double f9a = 4.0*pcnst_12*f7a;
      const double f10 = pcnst_32;

      const double  c_K1   =  1;
      const double  c_Na   =  1 *c20;
      const double  c_bNa  =  1*c21;
      const double  c_CaL  =  1;
      const double  c_bCa  =  1*c7;
      const double  c_NaCa =  1;
      const double  c_pCa  =  1*c24;
      const double  c_pK   =  1;
      const double  c_up   = 1*c44;
      const double  c_leak = 1*c43;
      const double  c_xfer = 1*c23;
      const double  c_rel  = 1*c40;

      const double c_NaK = P_NaK;
      const double c_Ks  = g_Ks;
      const double c_Kr  = g_Kr;
      const double c_to  = g_to;
      const double c_NaL = g_NaL;

      //double *states = cell[ii].state;
      double dVR = 0.0;
      //double itmp0,itmp5,itmp6;
      double I_K1,I_Kr,I_Ks,I_Na,I_bNa,I_CaL,I_bCa,I_to,I_NaK,I_NaCa,I_pCa,I_pK,I_NaL,I_leak,I_up,I_rel,I_xfer;
      double I_sum,I_delta;

      //  Update Ca concentration;
      {
         double fv1; {
            const double x = thisVm;
            const int l=1,m=6;
            const double a[]={1.5757133041269722e-04,2.4132516902459823e-06,1.4605316459098998e-08,6.0427222047915772e-11,2.6019285647913959e-13,6.0217000704020008e-16,1.0000000000000000e+00
            };
            double sum1 = 0;
            for (int j = m-1; j >= 0; j--)
               sum1 = a[j] + x*sum1;

            fv1 = sum1;
         }

         double fv2; {
            const double x = thisVm;
            const int l=1,m=8;
            const double a[]={5.4050480762039217e+02,-1.1953189161096423e+01,1.1865349408070526e-01,-6.0092610681801851e-04,4.0077214280113720e-07,2.9421099617665475e-08,-2.5528046472894199e-10,-1.5237049656618966e-12,1.0000000000000000e+00
            };
            double sum1 = 0;
            for (int j = m-1; j >= 0; j--)
               sum1 = a[j] + x*sum1;

            fv2 = sum1;
         }
         double x1 = _Ca_i * c26;
         double x2 = SQ(_Ca_i * c27);
         double x3 = SQ(_Ca_i * c28 + c29);
         double sigm1 = sigm(x1);
         double sigm2 = sigm(x2);
         double sigm3 = sigm(x3);
         double dV3 = thisVm - 0.5 * c3 * log(_Ca_i) - c8;

         I_NaCa = c_NaCa * (CUBE(_Na_i) * fv1 - _Ca_i * fv2);

         I_bCa = c_bCa * dV3;
         I_pCa = c_pCa * sigm1;
         I_up  =  c44 * sigm2;
         I_leak=  c43 * (_Ca_SR - _Ca_i);
         I_xfer = c23 * (_Ca_ss - _Ca_i);
         I_delta = I_leak - I_up; // I_detal = -itmp5
         I_sum =   I_bCa + I_pCa; // itmp4 =I_bCa+I_pCa;
         stateData[Ca_i_off*nCells_+ii]   = _Ca_i + (dt * c9) * (sigm3 * (0.5 * I_sum - I_NaCa + I_xfer * c15 + I_delta*c16));
         dVR  -= I_sum;
         //if (ii %4 ==0) printf("\n%d dVR=%14.12f ", ii,dVR);
         //else printf("%14.12f ", ii,I_sum);
      }

      double iK;
      //  Update K and Na concentrations;
      {
         double _K_i  = c9 * (_dVK_i - thisVm);
         double x0 = _Na_i * c25;

         double dV0 = thisVm - c3*log(_K_i) - c5;
         double dV1 = thisVm - c3*log(_Na_i) - c4;
         double dV2 = dV0 - c3 * logSeries(c2 * _Na_i/_K_i) + c5 - c6; // Assumption:  c2*_Na_i/_K_i is small;

         double fv0; {
            const double x = thisVm;
            const int l=1,m=10;
            const double a[]={-1.4549948900403902e+00,-2.2427240624934000e-03,2.8701981811073087e-05,-2.9352599282738483e-07,1.8687813196724116e-09,-6.1533240726750350e-12,-3.0715514390797755e-14,2.1068763323371155e-15,-1.5741479441621198e-17,-1.6908987107052300e-19, 1.0000000000000000e+00
            };
            double sum1 = 0;
            for (int j = m-1; j >= 0; j--)
               sum1 = a[j] + x*sum1;

            fv0 = sum1;

         }
         double fv5; {
            const double x = thisVm;
            const int l=5,m=8;
            const double a[]={2.1983170611730303e-04,2.1698309464995643e-05,9.8317099143409939e-07,2.6134710139327942e-08,4.3109819211796886e-10,4.3229574960011243e-12,2.4019472580987271e-14,5.6472223920064411e-17,1.0000000000000000e+00,-6.6157324225412789e-02,2.0076473883411361e-03,-2.6510671504853674e-05,2.2430697618084790e-07
            };
            double sum1 = 0;
            for (int j = m-1; j >= 0; j--)
               sum1 = a[j] + x*sum1;

            double sum2 = 0;
            int k = m + l - 1;
            for (int j = k; j >= m; j--)
               sum2 = a[j] + x * sum2;
            fv5 = sum1/sum2;
         }

         double fv6; {
            const double x = dV0;
            const int l=11,m=7;
            const double a[]={2.9794200055170916e-01,-8.8983788808079303e-03,7.9933972784803960e-03,-1.7599632267422598e-04,1.5074733144798214e-06,-5.9120168317614184e-09,8.9356468869962987e-12,1.0000000000000000e+00,2.1722631714279653e-01,2.1762877423819233e-02,1.5923509782789353e-03,-7.4897050385073499e-05,3.3873719083779318e-06,-5.3689489360972982e-08,8.8474419958442444e-10,-4.4602715802630215e-12,2.0785335146560804e-14,1.9580847130746169e-16
            };
            double sum1 = 0;
            for (int j = m-1; j >= 0; j--)
               sum1 = a[j] + x*sum1;

            double sum2 = 0;
            int k = m + l - 1;
            for (int j = k; j >= m; j--)
               sum2 = a[j] + x * sum2;
            fv6 = sum1/sum2;
         }

         I_NaK = c_NaK*sigm(x0) * fv0;                          // renamed itmpA to I_NaK
         I_pK  = c_pK * fv5 * dV0;
         I_K1  = c_K1 * fv6 * dV0;
         I_to  = c_to * rGate * sGate * dV0 ;
         I_Kr  = c_Kr * Xr1Gate * Xr2Gate * dV0;
         I_Na  = c_Na * CUBE(mGate) * hGate * jGate * dV1;
         I_NaL = c_NaL * CUBE(mGate) * jLGate * dV1;
         I_bNa = c_bNa * dV1;
         I_Ks  = c_Ks * SQ(XsGate) * dV2;
         
         double iNa =  3 * I_NaCa - 1.5 * I_NaK + I_Na + I_bNa + I_NaL;
         iK =  I_Ks + I_NaK + I_pK + I_K1 + I_to + I_Kr;


         stateData[dVK_i_off*nCells_+ii] += dt * iK;
         dVR    -=  iNa + iK - 2.0 * I_NaCa;
         stateData[Na_i_off*nCells_+ii] = _Na_i + (dt * c9) * iNa;
      }


      //  Update Ca_SS, Ca_SR, R_prime concentrations and fCass gate;
      {
         double fv3; {
            const double x = thisVm;
            const int l=1,m=11;
            const double a[]={1.0396540551771880e+00 ,4.6058010729117471e-02,7.9107591823998407e-04,4.3024156297446584e-06,-4.9702295945293119e-08,-6.9663904851401655e-10,2.4991578397134279e-12,7.8166616197552177e-14,1.3685444078468510e-16,-3.6584729114311930e-18,-1.7450648728130940e-20,1.0000000000000000e+00
            };
            double sum1 = 0;
            for (int j = m-1; j >= 0; j--)
               sum1 = a[j] + x*sum1;

            fv3 = sum1;
         }

         double fv4; {
            const double x = thisVm;
            const int l=1,m=9;
            const double a[]={2.5568410537089715e+01,-7.8144592191919016e-01,6.3170153846905984e-03,3.4110729724745464e-05,-3.7991861595738163e-07,-5.2236184692414528e-09,1.1857805262424718e-11,4.5506699371915196e-13,1.8913634488808377e-15,1.0000000000000000e+00
            };
            double sum1 = 0;
            for (int j = m-1; j >= 0; j--)
               sum1 = a[j] + x*sum1;

            fv4 = sum1;
         }
         double x4 = SQ(_Ca_SR * c30);
         double x5 = SQ(_Ca_SR * c31 + c32);
         double x6 = SQ(_Ca_ss * c33 + c34);
         double sigm4 = sigm(x4);
         double sigm5 = sigm(x5);
         double sigm6 = sigm(x6);
         double tmp8  = (c18 + c19 * sigm4); //Sigm4
         double tmp9  = tmp8 * _Ca_ss + c36;

         I_CaL = c_CaL * dGate * fGate * f2Gate * _fCass * (fv4 - _Ca_ss * fv3);  // renamed itmp1 to I_CaL

         double O = SQ(_Ca_ss) * _R_prime / (tmp8 * c17 + SQ(_Ca_ss));
         I_rel =c40 * O * (_Ca_SR - _Ca_ss);

         stateData[Ca_ss_off*nCells_+ii]   = _Ca_ss   + (dt * c9) * sigm6 * (I_xfer + I_rel * c14 + I_CaL * c13);
         stateData[Ca_sr_off*nCells_+ii]   = _Ca_SR   - (dt * c9) * sigm5 * (I_delta + I_rel);
         stateData[R_prime_off*nCells_+ii] = _R_prime + (dt * c9) * (c36 - tmp9 * _R_prime);

         //#if fCassForm == TT06
         double t1 = 1.0/(1.0 + SQ(20 * _Ca_ss));
         double mhu = 0.600000 * t1 + 0.4000000;
         double tauR = 1.0/(80.0 * t1 + 2.0);
         //#endif

         /*#if  fCassForm == RICE
         double xCa_ss = 1000*_Ca_ss;
         double xCa_i  = 1000*_Ca_i ;
         double mhu    = 0.6/(1.0 + xCa_ss * sqrt(xCa_ss) + xCa_i * sqrt(xCa_i)) + 0.4;
         double tauR   = 0.005/mhu;
         #endif
         */
         stateData[fCass_off*nCells_+ii]   = _fCass   + dt*(mhu - _fCass) * tauR;

         dVR += I_CaL;
      }


      dVm[ii] = dVR;
      
#define ratPolyGate()                           \
      double sum1,sum2;                         \
      const double x = Vm[ii];                  \
                                                \
      sum1 = 0;                                 \
      for (int j = Mhu_m-1; j >= 0; j--)        \
         sum1 = Mhu_a[j] + x*sum1;              \
                                                \
      sum2 = 0;                                 \
      int k = Mhu_m + Mhu_l - 1;                \
      for (int j = k; j >= Mhu_m; j--)          \
         sum2 = Mhu_a[j] + x * sum2;            \
      double mhu = sum1/sum2;                   \
                                                \
      sum1 = 0;                                 \
      for (int j = Tau_m-1; j >= 0; j--)        \
         sum1 = Tau_a[j] + x*sum1;              \
                                                \
      sum2 = 0;                                 \
      k = Tau_m + Tau_l - 1;                    \
      for (int j = k; j >= Tau_m; j--)          \
         sum2 = Tau_a[j] + x * sum2;            \
      double tauR = sum1/sum2*dt

#define RushLarsen(name) stateData[name##_off*nCells_+ii] -= (mhu - stateData[name##_off*nCells_+ii])*expm1(-tauR);
#define ForwardEuler(name) stateData[name##_off*nCells_+ii] += (mhu - stateData[name##_off*nCells_+ii])*tauR
      //0
      {
      const int Mhu_l = 10;
      const int Mhu_m = 5;
      const double Mhu_a[] = { 9.9632117206253790e-01,  4.0825738726469545e-02,  6.3401613233199589e-04,  4.4158436861700431e-06,  1.1622058324043520e-08,  1.0000000000000000e+00,  4.0568375699663400e-02,  6.4216825832642788e-04,  4.2661664422410096e-06,  1.3559930396321903e-08, -1.3573468728873069e-11, -4.2594802366702580e-13,  7.6779952208246166e-15,  1.4260675804433780e-16, -2.6656212072499249e-18};
      const int Tau_l = 1;
      const int Tau_m = 18;
      const double Tau_a[] = {1.7765862602413648e+01,  5.0010202770602419e-02, -7.8002064070783474e-04, -6.9399661775931530e-05,  1.6936588308244311e-06,  5.4629017090963798e-07, -1.3805420990037933e-08, -8.0678945216155694e-10,  1.6209833004622630e-11,  6.5130101230170358e-13, -6.9931705949674988e-15, -3.1161210504114690e-16,  5.0166191902609083e-19,  7.8608831661430381e-20,  4.3936315597226053e-22, -7.0535966258003289e-24, -9.0473475495087118e-26, -2.9878427692323621e-28,  1.0000000000000000e+00};

      ratPolyGate();
      RushLarsen(mGate);
      }
      //1
      {
         const int Mhu_l = 7;
         const int Mhu_m = 9;
         const double Mhu_a[] = {
            3.1263044718396443e-08,6.1534783640833612e-09,-3.7646458226792461e-10,-4.2451512350227289e-11,  1.5216916689168596e-12,  3.9843764909619204e-14, -7.7509942611798927e-16, -3.3646392579644030e-17,  5.7566741989713017e-19,  1.0000000000000000e+00,  6.4169311162380532e-02,  1.7346341071731517e-03,  2.5211977280197058e-05,  2.0770353677292379e-07,  9.2161835229419180e-10,  1.7368930279808849e-12
         };
         const int Tau_l = 11;
         const int Tau_m = 11;
         const double Tau_a[] = {
            4.2835355918249460e+00,  6.4588153679793869e-01,  4.4298327200183743e-02,  1.8236040756549200e-03,  5.0201035822685939e-05,  9.7279248118950919e-07,  1.3522528070830968e-08,  1.3356095236104801e-10,  8.9640565060229970e-13,  3.6781270948185535e-15,  6.9724701901279024e-18,  1.0000000000000000e+00,  1.2584512299220066e-01,  7.7047331771567910e-03,  3.0347083105033588e-04,  8.3931798507141038e-06,  1.6509899343992701e-07,  2.2953983633288216e-09,  2.2426481120489624e-11,  1.5049336925781968e-13,  6.2806501303798073e-16,  1.2085940888554432e-18
         };
         ratPolyGate();
         ForwardEuler(hGate);
      }
      //2
      {
         const int Mhu_l = 7;
         const int Mhu_m = 9;
         const double Mhu_a[] = {
            3.1263044718396443e-08,  6.1534783640833612e-09, -3.7646458226792461e-10, -4.2451512350227289e-11,  1.5216916689168596e-12,  3.9843764909619204e-14, -7.7509942611798927e-16, -3.3646392579644030e-17,  5.7566741989713017e-19,  1.0000000000000000e+00,  6.4169311162380532e-02,  1.7346341071731517e-03,  2.5211977280197058e-05,  2.0770353677292379e-07,  9.2161835229419180e-10,  1.7368930279808849e-12
         };
         const int Tau_l = 1;
         const int Tau_m = 13;
         const double Tau_a[] = {
            5.7685269377162662e-01,  3.5119791892517939e-02,  9.5403324816453923e-04,  1.8482121521375892e-05,  2.9493061828484874e-07,  2.6326436079978290e-09,  9.2707225823646994e-12,  3.6387910084447416e-13,  8.2748066828173235e-15,  2.3613021276064166e-17, -8.1293559809646280e-19, -7.9897103307330249e-21, -2.1875746179888087e-23,  1.0000000000000000e+00
         };
         ratPolyGate();
         ForwardEuler(jGate);
      }
      //3
      {
         const int Mhu_l = 5;
         const int Mhu_m = 8;
         const double Mhu_a[] = {
            9.7620687511302562e-01,  5.5687557065868629e-02,  1.3465134489445548e-03,  1.7240843920997653e-05,  1.1567866757841745e-07,  2.8292918822319700e-10, -8.2504917317015575e-13, -4.5498395933488446e-15,  1.0000000000000000e+00,  5.3646986850638727e-02,  1.4282138385673742e-03,  1.5252838752755205e-05,  1.4682559892184376e-07
         };
         const int Tau_l = 1;
         const int Tau_m = 13;
         const double Tau_a[] = {
            5.4609851558432825e-03,  4.3637206691412927e-04,  1.9034136985586222e-05,  5.4971744359155369e-07,  1.2019814013692350e-08,  2.0976070086137667e-10,  3.0204999531382605e-12,  3.6013108734580622e-14,  3.9727446952762997e-16,  4.3797461895814199e-18,  4.2418783268539835e-20,  2.6379981854029984e-22,  7.6831345074599439e-25,  1.0000000000000000e+00
         };
         ratPolyGate();
         ForwardEuler(Xr1Gate);

      }
      //4
      {
         const int Mhu_l = 1;
         const int Mhu_m = 10;
         const double Mhu_a[] = {
            2.4919224554816769e-02, -1.0132133908609619e-03,  2.0115370077599340e-05, -2.4787509921554170e-07,  1.9157062546455413e-09, -1.1856625513158431e-11,  4.3373559510489725e-14,  2.9410875655984132e-15, -2.8535346738989432e-17, -2.7319284078381429e-19,  1.0000000000000000e+00
         };
         const int Tau_l = 1;
         const int Tau_m = 10;
         const double Tau_a[] = {
            3.2797890377243832e-01,  1.1724387906430438e-06,  3.7188870926079577e-05, -3.3513027526302832e-09,  7.4502277507770199e-09,  1.2888652361968656e-12,  8.0734422444730268e-13,  8.1028416357975461e-16, -1.9648667839939542e-18, -3.5560477657149278e-19,  1.0000000000000000e+00
         };
         ratPolyGate();
         ForwardEuler(Xr2Gate);

      }
      //5
      {
         const int Mhu_l = 5;
         const int Mhu_m = 5;
         const double Mhu_a[] = {
            5.8834144047581471e-01,  1.8554918971499884e-02,  2.3769494810239468e-04,  1.4558400579866515e-06,  3.5608483256093079e-09,  1.0000000000000000e+00,  2.1263634327928394e-03,  5.2685461356465763e-04, -1.2603117542017038e-06,  1.5082411859988069e-08
         };
         const int Tau_l = 6;
         const int Tau_m = 9;
         const double Tau_a[] = {
            1.2782667970990670e-03, -5.3480687790940232e-05,  4.1463175377037539e-06, -7.9427267607107334e-08,  2.3181534592617566e-09, -1.8597152793072996e-11,  3.4109130605190490e-13,  2.0323839081334259e-15,  4.3130861032059193e-18,  1.0000000000000000e+00,  4.9865573111495263e-03,  1.0343375562869115e-03,  1.3748473468950345e-05,  3.3217061243887820e-07, -2.7950790321181472e-09
         };
         ratPolyGate();
         ForwardEuler(XsGate);
      }
      //6
      {
         const int Mhu_l = 5;
         const int Mhu_m = 8;
         const double Mhu_a[] = {
            3.4450952356981063e-02,  3.2569778153723060e-03,  1.4015137277741737e-04,  3.5177311829244476e-06,  5.4744904248961321e-08,  5.1978297648799841e-10,  2.7509401026136440e-12,  6.2017518287576868e-15,  1.0000000000000000e+00, -6.6423744780337163e-02,  2.2756546586311602e-03, -3.0510539067819270e-05,  3.1936830897476799e-07
         };
         const int Tau_l = 5;
         const int Tau_m = 7;
         const double Tau_a[] = {
            2.1250755813646321e-01,  5.4375121806696314e-03,  1.2855334968507940e-04,  1.8435064622728614e-06,  2.2201936942212559e-08,  1.6208374761956396e-10,  7.0475957646256506e-13,  1.0000000000000000e+00, -1.1298206006665951e-02,  1.9108836371246739e-05,  3.8822652158898901e-06,  2.4266538879422829e-08
         };
         ratPolyGate();
         ForwardEuler(rGate);
      }
      //7
      {
         const int Mhu_l = 5;
         const int Mhu_m = 7;
         const double Mhu_a[] = {
            7.4398777641576486e-01,  4.4326834369266523e-02,  1.1359075050132548e-03,  1.5872707499276087e-05,  1.2622142062210935e-07,  5.3559505029241141e-10,  9.3621878808126576e-13,  1.0000000000000000e+00,  2.5428070655441278e-02,  1.7701065218545982e-03,  3.7348639823389864e-06,  2.5522660570009651e-07
         };
         const int Tau_l = 7;
         const int Tau_m = 10;
         const double Tau_a[] = {
            1.5068318180921829e+00,  7.8993621885004611e-02,  6.9687982613764387e-03,  1.1019128271535092e-04,  2.3993512027253813e-06,  1.1356515133023819e-08,  3.1836750643625883e-10, -1.6718952430054101e-12,  2.5579267994953318e-14,  1.0819050682737564e-16,  1.0000000000000000e+00, -6.8101965713816826e-02,  3.7430651166330204e-03,  6.6379282851072469e-07,  9.3098141480209233e-07,  3.2614589730112466e-08,  4.4763500756297451e-10
         };
         ratPolyGate();
         ForwardEuler(dGate);
      }
      //8
      {
         const int Mhu_l = 5;
         const int Mhu_m = 8;
         const double Mhu_a[] = {
            5.4304780098652715e-02, -4.4962960006268609e-03,  1.7191393153273102e-04, -3.9003819270270578e-06,  5.4530196576796994e-08, -4.1503002319748868e-10,  8.3481595187685469e-13,  5.8904944849510583e-15,  1.0000000000000000e+00,  5.2318145763546209e-02,  1.6278824780069224e-03,  1.6010969045398731e-05,  1.9936720027259198e-07
         };
         const int Tau_l = 13;
         const int Tau_m = 11;
         const double Tau_a[] = {
            8.7301971950165709e-03,  1.0619909606051185e-03,  7.4439522770957243e-05,  3.5189522703546544e-06,  1.2763275154107722e-07,  3.9329666177349803e-09,  1.0963505894631244e-10,  2.5405559223755679e-12,  4.1926387533661996e-14,  4.1426318225610267e-16,  1.8370569627175333e-18,  1.0000000000000000e+00,  5.3496360858059103e-02,  1.0576677271235941e-02,  6.0292981147301659e-04,  2.4470545201049391e-05,  8.8342556285406000e-07,  2.5756751722220191e-08,  5.6840520101297543e-10,  8.9563864907355775e-12,  8.8109045583454404e-14,  4.1790988048136637e-16,  3.0481930918230229e-19,  7.4195798494855779e-22
         };
         ratPolyGate();
         ForwardEuler(fGate);
      }
      //9
      {
         const int Mhu_l = 8;
         const int Mhu_m = 5;
         const double Mhu_a[] = {
            3.3448475117101473e-01,  1.7908684328425843e-02,  4.6473146231373311e-04,  4.9947394623057196e-06,  4.6013519663601120e-08,  1.0000000000000000e+00,  5.5442666103601546e-02,  1.3608399868819414e-03,  1.6306821023781025e-05,  1.2063672477786111e-07,  1.7672562299303354e-10, -5.4528541905004782e-13, -2.8928059365721159e-15
         };
         const int Tau_l = 11;
         const int Tau_m = 12;
         const double Tau_a[] = {
            3.0215260344315119e-02,  2.9905565403315244e-03,  1.8698675038604531e-04,  7.8302945745729333e-06,  2.4353851952274140e-07,  5.9716586064763031e-09,  1.4320157611618665e-10,  3.5798223309278773e-12,  7.3139244140308910e-14,  9.1797821675112518e-16,  5.5127688299509001e-18,  6.0160642283988931e-21,  1.0000000000000000e+00, -8.8631007579573415e-02,  5.6222208830099462e-03,  9.0331117197141856e-05,  8.5481830250918979e-07,  2.2780719183249479e-07,  7.0378379382184305e-09,  1.0211239994063565e-10,  1.8042686984162349e-12,  2.9931381582699879e-14,  2.1128053668752277e-16
         };
         ratPolyGate();
         ForwardEuler(f2Gate);
      }
      //10
      {
         const int Mhu_l = 6;
         const int Mhu_m = 12;
         const double Mhu_a[] = {
            3.6932725836170428e-09,  2.5736169972850485e-10, -7.5177270851662063e-11, -1.7783383158540549e-12,  2.3780959584764758e-13,  3.9995914648833690e-15, -2.6002507057823451e-16, -4.0173750386716292e-18,  1.1567747025479995e-19,  1.8190095403795220e-21, -1.7772014587013272e-23, -3.0513026345262148e-25,  1.0000000000000000e+00,  4.7761481528554264e-02,  9.1744037118706855e-04,  8.8560239615590481e-06,  4.2943755440989806e-08,  8.3656852902593096e-11
         };
         const int Tau_l = 1;
         const int Tau_m = 1;
         const double Tau_a[] = {
            1.4925373134328358e-03,  1.0000000000000000e+00
         };
         ratPolyGate();
         ForwardEuler(jLGate);
      }
      //11
      {
         const int Mhu_l = 7;
         const int Mhu_m = 8;
         const double Mhu_a[] = {
            3.6866939888900814e-03, -4.5574539113452475e-04,  2.6788146214637530e-05, -9.9102133439606097e-07,  2.5613910193807883e-08, -4.6400777680998212e-10,  5.2352368345974410e-12, -2.6731132916017718e-14,  1.0000000000000000e+00,  7.5601232522433295e-02,  2.5744765696710309e-03,  4.7389711229137498e-05,  5.8624244628409614e-07,  3.4873125981144344e-09,  2.0871941515527373e-11
         };
         const int Tau_l = 6;
         const int Tau_m = 7;
         const double Tau_a[] = {
            5.2008879668973697e-02,  2.9402318950257256e-03,  9.6193142778593343e-05,  1.9694354926101320e-06,  2.5766968301605511e-08,  1.9309877695248136e-10,  6.5471492087484267e-13,  1.0000000000000000e+00, -2.1715782138805460e-02,  2.0915626226951107e-03,  1.8073629807905746e-07,  2.1067307291025196e-07,  2.7434857789930993e-09
         };
         ratPolyGate();
         ForwardEuler(sGate);
      }   
}

void actualCalc(const double dt, const int nCells_, const double Vm[], const double iStim[], double dVm[], double* stateData) {
   const int blockSize=512;
   
   actualCalcCuda<<<dim3((nCells_+blockSize-1)/blockSize,1,1),dim3(blockSize,1,1)>>>(dt,nCells_,Vm,iStim,dVm,stateData);
}

const char* SimpleTT06_ker = "SimpleTT06_ker\n"
"\n"
"   \n"
"enum StateOffset {\n"
"   f2Gate_off,\n"
"   fGate_off,\n"
"   dGate_off,\n"
"   mGate_off,\n"
"   jGate_off,\n"
"   hGate_off,\n"
"   rGate_off,\n"
"   sGate_off,\n"
"   Xr1Gate_off,\n"
"   Xr2Gate_off,\n"
"   XsGate_off,\n"
"   jLGate_off,\n"
"   Na_i_off,\n"
"   Ca_i_off,\n"
"   Ca_ss_off,\n"
"   Ca_sr_off,\n"
"   fCass_off,\n"
"   dVK_i_off,\n"
"   R_prime_off,\n"
"   NUMSTATES\n"
"};\n"
"\n"
"#define CUBE(x) ((x)*(x)*(x))\n"
"#define SQ(x) ((x)*(x))\n"
"#define sigm(x)   ((x)/(1+(x)))\n"
"#define logSeries(x)    (log(1+(x)) )\n"
"\n"
"__global__ void actualCalcCudaJitify(const double dt, const int nCells_, const double Vm[], const double iStim[], double dVm[], double* stateData)\n"
"{   \n"
"   int ii = blockIdx.x*blockDim.x + threadIdx.x;\n"
"   if (ii > nCells_) { return; }\n"
"   //printf(\"CPU state = %p\\n\", &state);\n"
"   //printf(\"CPU state->Ca_i = %p\\n\", state->Ca_i);\n"
"   //printf(\"CPU state->nCells = %d\\n\", state->nCells);\n"
"   //printf(\"CPU istim[0] = %g\\n\", iStim[0]);\n"
"\n"
"      //printf(\"CPU state = %p\\n\", &state);\n"
"      //printf(\"CPU state->Ca_i = %p\\n\", state->Ca_i);\n"
"      //printf(\"CPU state->nCells = %d\\n\", state->nCells);\n"
"      //printf(\"GPU istim[0] = %g\\n\", iStim[0]);\n"
"\n"
"\n"
"      //set Vm\n"
"      const double thisVm = Vm[ii];\n"
"      const double istim = iStim[ii];\n"
"\n"
"      //set all state variables\n"
"      //EDIT_STATE\n"
"      const double f2Gate=stateData[f2Gate_off*nCells_+ii];\n"
"      const double fGate=stateData[fGate_off*nCells_+ii];\n"
"      const double dGate=stateData[dGate_off*nCells_+ii];\n"
"      const double mGate=stateData[mGate_off*nCells_+ii];\n"
"      const double jGate=stateData[jGate_off*nCells_+ii];\n"
"      const double hGate=stateData[hGate_off*nCells_+ii];\n"
"      const double rGate=stateData[rGate_off*nCells_+ii];\n"
"      const double sGate=stateData[sGate_off*nCells_+ii];\n"
"      const double Xr1Gate=stateData[Xr1Gate_off*nCells_+ii];\n"
"      const double Xr2Gate=stateData[Xr2Gate_off*nCells_+ii];\n"
"      const double XsGate=stateData[XsGate_off*nCells_+ii];\n"
"      const double jLGate=stateData[jLGate_off*nCells_+ii];\n"
"      const double _Na_i=stateData[Na_i_off*nCells_+ii];\n"
"      const double _Ca_i=stateData[Ca_i_off*nCells_+ii];\n"
"      const double _Ca_ss=stateData[Ca_ss_off*nCells_+ii];\n"
"      const double _Ca_SR=stateData[Ca_sr_off*nCells_+ii];\n"
"      const double _fCass=stateData[fCass_off*nCells_+ii];\n"
"      const double _dVK_i=stateData[dVK_i_off*nCells_+ii];\n"
"      const double _R_prime=stateData[R_prime_off*nCells_+ii];\n"
"\n"
"      //set per-cell flags\n"
"      //EDIT_PERCELL_FLAGS\n"
"      \n"
"      //set per-cell parameters\n"
"      //EDIT_PERCELL_PARAMETERS\n"
"\n"
"      const double pcnst_0 = 8314.472;\n"
"      const double pcnst_1 = 310;\n"
"      const double pcnst_2 = 96485.3415;\n"
"      const double pcnst_9 = 0.03;\n"
"      const double pcnst_10 = 5.4;\n"
"      const double pcnst_11 = 140;\n"
"      const double pcnst_12 = 2;\n"
"      const double c1 = pcnst_2/(pcnst_0*pcnst_1);\n"
"      const double c2 = pcnst_9;\n"
"      const double c3 = -1/c1;\n"
"\n"
"      //const double c4 = -c3*log(pcnst_11);\n"
"      //const double c5 = -c3*log(pcnst_10);\n"
"      //const double c6 = -c3*log(pcnst_10+pcnst_9*pcnst_11);\n"
"      //const double c8 = -0.5*c3*log(pcnst_12);\n"
"      //if (0) {\n"
"      //  printf(\"%.17g %.17g %.17g %.17g\\n\", c4,c5,c6,c8);\n"
"      //}\n"
"\n"
"      const double c4 = 132.00985294338352;\n"
"      const double c5 = 45.050058022436289;\n"
"      const double c6 = 60.420198081560486;\n"
"      const double c8 = 9.2582839417106104;\n"
"      const double P_NaK=2.724;\n"
"      const double g_Ks = 0.392;\n"
"      const double g_Kr = 0.153;\n"
"      const double g_to = 0.073;\n"
"      const double g_NaL=0.0;\n"
"\n"
"      const double pcnst_3 = 0.185;\n"
"      const double pcnst_4 = 0.016404;\n"
"      const double pcnst_5 = 10;\n"
"      const double pcnst_6 = 1000;\n"
"      const double pcnst_7 = 1;\n"
"      const double pcnst_8 = 52;\n"
"      const double pcnst_13 = 5.405;\n"
"      const double pcnst_14 = 0.153;\n"
"      const double pcnst_16 = 14.838;\n"
"      const double pcnst_17 = 0.00029;\n"
"      const double pcnst_18 = 0.0000398;\n"
"      const double pcnst_19 = 0.000592;\n"
"      const double pcnst_21 = 2.724;\n"
"      const double pcnst_22 = 1;\n"
"      const double pcnst_23 = 40;\n"
"      const double pcnst_24 = 1000;\n"
"      const double pcnst_25 = 0.1;\n"
"      const double pcnst_26 = 2.5;\n"
"      const double pcnst_27 = 0.35;\n"
"      const double pcnst_28 = 1.38;\n"
"      const double pcnst_29 = 87.5;\n"
"      const double pcnst_30 = 0.1238;\n"
"      const double pcnst_31 = 0.0005;\n"
"      const double pcnst_32 = 0.0146;\n"
"      const double pcnst_33 = 0.15;\n"
"      const double pcnst_34 = 0.045;\n"
"      const double pcnst_35 = 0.06;\n"
"      const double pcnst_36 = 0.005;\n"
"      const double pcnst_37 = 1.5;\n"
"      const double pcnst_38 = 2.5;\n"
"      const double pcnst_39 = 1;\n"
"      const double pcnst_40 = 0.102;\n"
"      const double pcnst_41 = 0.0038;\n"
"      const double pcnst_42 = 0.00025;\n"
"      const double pcnst_43 = 0.00036;\n"
"      const double pcnst_44 = 0.006375;\n"
"      const double pcnst_45 = 0.2;\n"
"      const double pcnst_46 = 0.001;\n"
"      const double pcnst_47 = 10;\n"
"      const double pcnst_48 = 0.3;\n"
"      const double pcnst_49 = 0.4;\n"
"      const double pcnst_50 = 0.00025;\n"
"      const double pcnst_51 = 0.001094;\n"
"      const double pcnst_52 = 0.00005468;\n"
"\n"
"      const double c9 = -pcnst_3/(pcnst_4*pcnst_2);\n"
"\n"
"\n"
"      const double c7 =  pcnst_19;\n"
"      //c11= pcsnt[14]*sqrt(pcnst_10/5.4);\n"
"      //const double c11= sqrt(pcnst_10/5.4);\n"
"      const double c11 = 1;\n"
"      //const double c12= pcnst_13*sqrt(pcnst_10/5.4);\n"
"      const double c12 = 5.4050000000000002;\n"
"      const double c13= pcnst_3/(2.0*pcnst_52*pcnst_2*c9);\n"
"      const double c14 = pcnst_51/pcnst_52;\n"
"      const double c15 = -pcnst_52/pcnst_4;\n"
"      const double c16 = pcnst_51/pcnst_4;\n"
"      const double c17 = pcnst_35/(pcnst_33*pcnst_34)*c9;\n"
"      const double c18 = pcnst_34*pcnst_38/c9 ;\n"
"      const double c19  = -pcnst_34*(pcnst_38-pcnst_39)/c9;\n"
"      const double c20  = pcnst_16;\n"
"      const double c21  = pcnst_17;\n"
"      const double c22  = 1/c9;\n"
"      const double c23  = pcnst_41/(c15*c9);\n"
"      const double c24  =  pcnst_30;\n"
"      const double c25  =  1.0/pcnst_23;\n"
"      const double c26  =  1.0/pcnst_31;\n"
"      const double c27  =  1.0/pcnst_42;\n"
"      //const double c28  =  1.0/sqrt(pcnst_45*pcnst_46);\n"
"      const double c28 = 70.710678118654755;\n"
"      const double c29  =  pcnst_46*c28;\n"
"      const double c30  =  1.0/pcnst_37;\n"
"      //const double c31  =  1.0/sqrt(pcnst_47*pcnst_48);\n"
"      const double c31 = 0.57735026918962584;\n"
"      const double c32  =  pcnst_48*c31;\n"
"      //const double c33  =  1.0/sqrt(pcnst_49*pcnst_50);\n"
"      const double c33 = 100;\n"
"      const double c34  =  pcnst_50*c33;\n"
"      const double c36  =  pcnst_36/c9;\n"
"      const double c40  =  pcnst_40/c9;\n"
"      const double c43  =  pcnst_43/c9;\n"
"      const double c44  =  pcnst_44/c9;\n"
"\n"
"      const double f1 = c1;\n"
"      const double f2 =  -2.0*pcnst_10/(pcnst_10+pcnst_22);\n"
"      const double f3 =  ((CUBE(pcnst_29)+CUBE(pcnst_11))*(pcnst_28+pcnst_12))/(pcnst_24*pcnst_12);\n"
"      const double f4 =  f3*pcnst_25;\n"
"   \n"
"      const double f5 =  pcnst_27*f1;\n"
"      const double f6 =  (CUBE(pcnst_11)*pcnst_26/pcnst_12);\n"
"      const double f7 = pcnst_18*pcnst_2*f1;\n"
"      const double f7a = 0.5*pcnst_18*pcnst_2;\n"
"      const double f9  = 4.0*pcnst_12;\n"
"      const double f9a = 4.0*pcnst_12*f7a;\n"
"      const double f10 = pcnst_32;\n"
"\n"
"      const double  c_K1   =  1;\n"
"      const double  c_Na   =  1 *c20;\n"
"      const double  c_bNa  =  1*c21;\n"
"      const double  c_CaL  =  1;\n"
"      const double  c_bCa  =  1*c7;\n"
"      const double  c_NaCa =  1;\n"
"      const double  c_pCa  =  1*c24;\n"
"      const double  c_pK   =  1;\n"
"      const double  c_up   = 1*c44;\n"
"      const double  c_leak = 1*c43;\n"
"      const double  c_xfer = 1*c23;\n"
"      const double  c_rel  = 1*c40;\n"
"\n"
"      const double c_NaK = P_NaK;\n"
"      const double c_Ks  = g_Ks;\n"
"      const double c_Kr  = g_Kr;\n"
"      const double c_to  = g_to;\n"
"      const double c_NaL = g_NaL;\n"
"\n"
"      //double *states = cell[ii].state;\n"
"      double dVR = 0.0;\n"
"      //double itmp0,itmp5,itmp6;\n"
"      double I_K1,I_Kr,I_Ks,I_Na,I_bNa,I_CaL,I_bCa,I_to,I_NaK,I_NaCa,I_pCa,I_pK,I_NaL,I_leak,I_up,I_rel,I_xfer;\n"
"      double I_sum,I_delta;\n"
"\n"
"      //  Update Ca concentration;\n"
"      {\n"
"         double fv1; {\n"
"            const double x = thisVm;\n"
"            const int l=1,m=6;\n"
"            const double a[]={1.5757133041269722e-04,2.4132516902459823e-06,1.4605316459098998e-08,6.0427222047915772e-11,2.6019285647913959e-13,6.0217000704020008e-16,1.0000000000000000e+00\n"
"            };\n"
"            double sum1 = 0;\n"
"            for (int j = m-1; j >= 0; j--)\n"
"               sum1 = a[j] + x*sum1;\n"
"\n"
"            fv1 = sum1;\n"
"         }\n"
"\n"
"         double fv2; {\n"
"            const double x = thisVm;\n"
"            const int l=1,m=8;\n"
"            const double a[]={5.4050480762039217e+02,-1.1953189161096423e+01,1.1865349408070526e-01,-6.0092610681801851e-04,4.0077214280113720e-07,2.9421099617665475e-08,-2.5528046472894199e-10,-1.5237049656618966e-12,1.0000000000000000e+00\n"
"            };\n"
"            double sum1 = 0;\n"
"            for (int j = m-1; j >= 0; j--)\n"
"               sum1 = a[j] + x*sum1;\n"
"\n"
"            fv2 = sum1;\n"
"         }\n"
"         double x1 = _Ca_i * c26;\n"
"         double x2 = SQ(_Ca_i * c27);\n"
"         double x3 = SQ(_Ca_i * c28 + c29);\n"
"         double sigm1 = sigm(x1);\n"
"         double sigm2 = sigm(x2);\n"
"         double sigm3 = sigm(x3);\n"
"         double dV3 = thisVm - 0.5 * c3 * log(_Ca_i) - c8;\n"
"\n"
"         I_NaCa = c_NaCa * (CUBE(_Na_i) * fv1 - _Ca_i * fv2);\n"
"\n"
"         I_bCa = c_bCa * dV3;\n"
"         I_pCa = c_pCa * sigm1;\n"
"         I_up  =  c44 * sigm2;\n"
"         I_leak=  c43 * (_Ca_SR - _Ca_i);\n"
"         I_xfer = c23 * (_Ca_ss - _Ca_i);\n"
"         I_delta = I_leak - I_up; // I_detal = -itmp5\n"
"         I_sum =   I_bCa + I_pCa; // itmp4 =I_bCa+I_pCa;\n"
"         stateData[Ca_i_off*nCells_+ii]   = _Ca_i + (dt * c9) * (sigm3 * (0.5 * I_sum - I_NaCa + I_xfer * c15 + I_delta*c16));\n"
"         dVR  -= I_sum;\n"
"         //if (ii %4 ==0) printf(\"\\n%d dVR=%14.12f \", ii,dVR);\n"
"         //else printf(\"%14.12f \", ii,I_sum);\n"
"      }\n"
"\n"
"      double iK;\n"
"      //  Update K and Na concentrations;\n"
"      {\n"
"         double _K_i  = c9 * (_dVK_i - thisVm);\n"
"         double x0 = _Na_i * c25;\n"
"\n"
"         double dV0 = thisVm - c3*log(_K_i) - c5;\n"
"         double dV1 = thisVm - c3*log(_Na_i) - c4;\n"
"         double dV2 = dV0 - c3 * logSeries(c2 * _Na_i/_K_i) + c5 - c6; // Assumption:  c2*_Na_i/_K_i is small;\n"
"\n"
"         double fv0; {\n"
"            const double x = thisVm;\n"
"            const int l=1,m=10;\n"
"            const double a[]={-1.4549948900403902e+00,-2.2427240624934000e-03,2.8701981811073087e-05,-2.9352599282738483e-07,1.8687813196724116e-09,-6.1533240726750350e-12,-3.0715514390797755e-14,2.1068763323371155e-15,-1.5741479441621198e-17,-1.6908987107052300e-19, 1.0000000000000000e+00\n"
"            };\n"
"            double sum1 = 0;\n"
"            for (int j = m-1; j >= 0; j--)\n"
"               sum1 = a[j] + x*sum1;\n"
"\n"
"            fv0 = sum1;\n"
"\n"
"         }\n"
"         double fv5; {\n"
"            const double x = thisVm;\n"
"            const int l=5,m=8;\n"
"            const double a[]={2.1983170611730303e-04,2.1698309464995643e-05,9.8317099143409939e-07,2.6134710139327942e-08,4.3109819211796886e-10,4.3229574960011243e-12,2.4019472580987271e-14,5.6472223920064411e-17,1.0000000000000000e+00,-6.6157324225412789e-02,2.0076473883411361e-03,-2.6510671504853674e-05,2.2430697618084790e-07\n"
"            };\n"
"            double sum1 = 0;\n"
"            for (int j = m-1; j >= 0; j--)\n"
"               sum1 = a[j] + x*sum1;\n"
"\n"
"            double sum2 = 0;\n"
"            int k = m + l - 1;\n"
"            for (int j = k; j >= m; j--)\n"
"               sum2 = a[j] + x * sum2;\n"
"            fv5 = sum1/sum2;\n"
"         }\n"
"\n"
"         double fv6; {\n"
"            const double x = dV0;\n"
"            const int l=11,m=7;\n"
"            const double a[]={2.9794200055170916e-01,-8.8983788808079303e-03,7.9933972784803960e-03,-1.7599632267422598e-04,1.5074733144798214e-06,-5.9120168317614184e-09,8.9356468869962987e-12,1.0000000000000000e+00,2.1722631714279653e-01,2.1762877423819233e-02,1.5923509782789353e-03,-7.4897050385073499e-05,3.3873719083779318e-06,-5.3689489360972982e-08,8.8474419958442444e-10,-4.4602715802630215e-12,2.0785335146560804e-14,1.9580847130746169e-16\n"
"            };\n"
"            double sum1 = 0;\n"
"            for (int j = m-1; j >= 0; j--)\n"
"               sum1 = a[j] + x*sum1;\n"
"\n"
"            double sum2 = 0;\n"
"            int k = m + l - 1;\n"
"            for (int j = k; j >= m; j--)\n"
"               sum2 = a[j] + x * sum2;\n"
"            fv6 = sum1/sum2;\n"
"         }\n"
"\n"
"         I_NaK = c_NaK*sigm(x0) * fv0;                          // renamed itmpA to I_NaK\n"
"         I_pK  = c_pK * fv5 * dV0;\n"
"         I_K1  = c_K1 * fv6 * dV0;\n"
"         I_to  = c_to * rGate * sGate * dV0 ;\n"
"         I_Kr  = c_Kr * Xr1Gate * Xr2Gate * dV0;\n"
"         I_Na  = c_Na * CUBE(mGate) * hGate * jGate * dV1;\n"
"         I_NaL = c_NaL * CUBE(mGate) * jLGate * dV1;\n"
"         I_bNa = c_bNa * dV1;\n"
"         I_Ks  = c_Ks * SQ(XsGate) * dV2;\n"
"         \n"
"         double iNa =  3 * I_NaCa - 1.5 * I_NaK + I_Na + I_bNa + I_NaL;\n"
"         iK =  I_Ks + I_NaK + I_pK + I_K1 + I_to + I_Kr;\n"
"\n"
"\n"
"         stateData[dVK_i_off*nCells_+ii] += dt * iK;\n"
"         dVR    -=  iNa + iK - 2.0 * I_NaCa;\n"
"         stateData[Na_i_off*nCells_+ii] = _Na_i + (dt * c9) * iNa;\n"
"      }\n"
"\n"
"\n"
"      //  Update Ca_SS, Ca_SR, R_prime concentrations and fCass gate;\n"
"      {\n"
"         double fv3; {\n"
"            const double x = thisVm;\n"
"            const int l=1,m=11;\n"
"            const double a[]={1.0396540551771880e+00 ,4.6058010729117471e-02,7.9107591823998407e-04,4.3024156297446584e-06,-4.9702295945293119e-08,-6.9663904851401655e-10,2.4991578397134279e-12,7.8166616197552177e-14,1.3685444078468510e-16,-3.6584729114311930e-18,-1.7450648728130940e-20,1.0000000000000000e+00\n"
"            };\n"
"            double sum1 = 0;\n"
"            for (int j = m-1; j >= 0; j--)\n"
"               sum1 = a[j] + x*sum1;\n"
"\n"
"            fv3 = sum1;\n"
"         }\n"
"\n"
"         double fv4; {\n"
"            const double x = thisVm;\n"
"            const int l=1,m=9;\n"
"            const double a[]={2.5568410537089715e+01,-7.8144592191919016e-01,6.3170153846905984e-03,3.4110729724745464e-05,-3.7991861595738163e-07,-5.2236184692414528e-09,1.1857805262424718e-11,4.5506699371915196e-13,1.8913634488808377e-15,1.0000000000000000e+00\n"
"            };\n"
"            double sum1 = 0;\n"
"            for (int j = m-1; j >= 0; j--)\n"
"               sum1 = a[j] + x*sum1;\n"
"\n"
"            fv4 = sum1;\n"
"         }\n"
"         double x4 = SQ(_Ca_SR * c30);\n"
"         double x5 = SQ(_Ca_SR * c31 + c32);\n"
"         double x6 = SQ(_Ca_ss * c33 + c34);\n"
"         double sigm4 = sigm(x4);\n"
"         double sigm5 = sigm(x5);\n"
"         double sigm6 = sigm(x6);\n"
"         double tmp8  = (c18 + c19 * sigm4); //Sigm4\n"
"         double tmp9  = tmp8 * _Ca_ss + c36;\n"
"\n"
"         I_CaL = c_CaL * dGate * fGate * f2Gate * _fCass * (fv4 - _Ca_ss * fv3);  // renamed itmp1 to I_CaL\n"
"\n"
"         double O = SQ(_Ca_ss) * _R_prime / (tmp8 * c17 + SQ(_Ca_ss));\n"
"         I_rel =c40 * O * (_Ca_SR - _Ca_ss);\n"
"\n"
"         stateData[Ca_ss_off*nCells_+ii]   = _Ca_ss   + (dt * c9) * sigm6 * (I_xfer + I_rel * c14 + I_CaL * c13);\n"
"         stateData[Ca_sr_off*nCells_+ii]   = _Ca_SR   - (dt * c9) * sigm5 * (I_delta + I_rel);\n"
"         stateData[R_prime_off*nCells_+ii] = _R_prime + (dt * c9) * (c36 - tmp9 * _R_prime);\n"
"\n"
"         //#if fCassForm == TT06\n"
"         double t1 = 1.0/(1.0 + SQ(20 * _Ca_ss));\n"
"         double mhu = 0.600000 * t1 + 0.4000000;\n"
"         double tauR = 1.0/(80.0 * t1 + 2.0);\n"
"         //#endif\n"
"\n"
"         /*#if  fCassForm == RICE\n"
"         double xCa_ss = 1000*_Ca_ss;\n"
"         double xCa_i  = 1000*_Ca_i ;\n"
"         double mhu    = 0.6/(1.0 + xCa_ss * sqrt(xCa_ss) + xCa_i * sqrt(xCa_i)) + 0.4;\n"
"         double tauR   = 0.005/mhu;\n"
"         #endif\n"
"         */\n"
"         stateData[fCass_off*nCells_+ii]   = _fCass   + dt*(mhu - _fCass) * tauR;\n"
"\n"
"         dVR += I_CaL;\n"
"      }\n"
"\n"
"\n"
"      dVm[ii] = dVR;\n"
"      \n"
"#define ratPolyGate()                           \\\n"
"      double sum1,sum2;                         \\\n"
"      const double x = Vm[ii];                  \\\n"
"                                                \\\n"
"      sum1 = 0;                                 \\\n"
"      for (int j = Mhu_m-1; j >= 0; j--)        \\\n"
"         sum1 = Mhu_a[j] + x*sum1;              \\\n"
"                                                \\\n"
"      sum2 = 0;                                 \\\n"
"      int k = Mhu_m + Mhu_l - 1;                \\\n"
"      for (int j = k; j >= Mhu_m; j--)          \\\n"
"         sum2 = Mhu_a[j] + x * sum2;            \\\n"
"      double mhu = sum1/sum2;                   \\\n"
"                                                \\\n"
"      sum1 = 0;                                 \\\n"
"      for (int j = Tau_m-1; j >= 0; j--)        \\\n"
"         sum1 = Tau_a[j] + x*sum1;              \\\n"
"                                                \\\n"
"      sum2 = 0;                                 \\\n"
"      k = Tau_m + Tau_l - 1;                    \\\n"
"      for (int j = k; j >= Tau_m; j--)          \\\n"
"         sum2 = Tau_a[j] + x * sum2;            \\\n"
"      double tauR = sum1/sum2*dt\n"
"\n"
"#define RushLarsen(name) stateData[name##_off*nCells_+ii] -= (mhu - stateData[name##_off*nCells_+ii])*expm1(-tauR);\n"
"#define ForwardEuler(name) stateData[name##_off*nCells_+ii] += (mhu - stateData[name##_off*nCells_+ii])*tauR\n"
"      //0\n"
"      {\n"
"      const int Mhu_l = 10;\n"
"      const int Mhu_m = 5;\n"
"      const double Mhu_a[] = { 9.9632117206253790e-01,  4.0825738726469545e-02,  6.3401613233199589e-04,  4.4158436861700431e-06,  1.1622058324043520e-08,  1.0000000000000000e+00,  4.0568375699663400e-02,  6.4216825832642788e-04,  4.2661664422410096e-06,  1.3559930396321903e-08, -1.3573468728873069e-11, -4.2594802366702580e-13,  7.6779952208246166e-15,  1.4260675804433780e-16, -2.6656212072499249e-18};\n"
"      const int Tau_l = 1;\n"
"      const int Tau_m = 18;\n"
"      const double Tau_a[] = {1.7765862602413648e+01,  5.0010202770602419e-02, -7.8002064070783474e-04, -6.9399661775931530e-05,  1.6936588308244311e-06,  5.4629017090963798e-07, -1.3805420990037933e-08, -8.0678945216155694e-10,  1.6209833004622630e-11,  6.5130101230170358e-13, -6.9931705949674988e-15, -3.1161210504114690e-16,  5.0166191902609083e-19,  7.8608831661430381e-20,  4.3936315597226053e-22, -7.0535966258003289e-24, -9.0473475495087118e-26, -2.9878427692323621e-28,  1.0000000000000000e+00};\n"
"\n"
"      ratPolyGate();\n"
"      RushLarsen(mGate);\n"
"      }\n"
"      //1\n"
"      {\n"
"         const int Mhu_l = 7;\n"
"         const int Mhu_m = 9;\n"
"         const double Mhu_a[] = {\n"
"            3.1263044718396443e-08,6.1534783640833612e-09,-3.7646458226792461e-10,-4.2451512350227289e-11,  1.5216916689168596e-12,  3.9843764909619204e-14, -7.7509942611798927e-16, -3.3646392579644030e-17,  5.7566741989713017e-19,  1.0000000000000000e+00,  6.4169311162380532e-02,  1.7346341071731517e-03,  2.5211977280197058e-05,  2.0770353677292379e-07,  9.2161835229419180e-10,  1.7368930279808849e-12\n"
"         };\n"
"         const int Tau_l = 11;\n"
"         const int Tau_m = 11;\n"
"         const double Tau_a[] = {\n"
"            4.2835355918249460e+00,  6.4588153679793869e-01,  4.4298327200183743e-02,  1.8236040756549200e-03,  5.0201035822685939e-05,  9.7279248118950919e-07,  1.3522528070830968e-08,  1.3356095236104801e-10,  8.9640565060229970e-13,  3.6781270948185535e-15,  6.9724701901279024e-18,  1.0000000000000000e+00,  1.2584512299220066e-01,  7.7047331771567910e-03,  3.0347083105033588e-04,  8.3931798507141038e-06,  1.6509899343992701e-07,  2.2953983633288216e-09,  2.2426481120489624e-11,  1.5049336925781968e-13,  6.2806501303798073e-16,  1.2085940888554432e-18\n"
"         };\n"
"         ratPolyGate();\n"
"         ForwardEuler(hGate);\n"
"      }\n"
"      //2\n"
"      {\n"
"         const int Mhu_l = 7;\n"
"         const int Mhu_m = 9;\n"
"         const double Mhu_a[] = {\n"
"            3.1263044718396443e-08,  6.1534783640833612e-09, -3.7646458226792461e-10, -4.2451512350227289e-11,  1.5216916689168596e-12,  3.9843764909619204e-14, -7.7509942611798927e-16, -3.3646392579644030e-17,  5.7566741989713017e-19,  1.0000000000000000e+00,  6.4169311162380532e-02,  1.7346341071731517e-03,  2.5211977280197058e-05,  2.0770353677292379e-07,  9.2161835229419180e-10,  1.7368930279808849e-12\n"
"         };\n"
"         const int Tau_l = 1;\n"
"         const int Tau_m = 13;\n"
"         const double Tau_a[] = {\n"
"            5.7685269377162662e-01,  3.5119791892517939e-02,  9.5403324816453923e-04,  1.8482121521375892e-05,  2.9493061828484874e-07,  2.6326436079978290e-09,  9.2707225823646994e-12,  3.6387910084447416e-13,  8.2748066828173235e-15,  2.3613021276064166e-17, -8.1293559809646280e-19, -7.9897103307330249e-21, -2.1875746179888087e-23,  1.0000000000000000e+00\n"
"         };\n"
"         ratPolyGate();\n"
"         ForwardEuler(jGate);\n"
"      }\n"
"      //3\n"
"      {\n"
"         const int Mhu_l = 5;\n"
"         const int Mhu_m = 8;\n"
"         const double Mhu_a[] = {\n"
"            9.7620687511302562e-01,  5.5687557065868629e-02,  1.3465134489445548e-03,  1.7240843920997653e-05,  1.1567866757841745e-07,  2.8292918822319700e-10, -8.2504917317015575e-13, -4.5498395933488446e-15,  1.0000000000000000e+00,  5.3646986850638727e-02,  1.4282138385673742e-03,  1.5252838752755205e-05,  1.4682559892184376e-07\n"
"         };\n"
"         const int Tau_l = 1;\n"
"         const int Tau_m = 13;\n"
"         const double Tau_a[] = {\n"
"            5.4609851558432825e-03,  4.3637206691412927e-04,  1.9034136985586222e-05,  5.4971744359155369e-07,  1.2019814013692350e-08,  2.0976070086137667e-10,  3.0204999531382605e-12,  3.6013108734580622e-14,  3.9727446952762997e-16,  4.3797461895814199e-18,  4.2418783268539835e-20,  2.6379981854029984e-22,  7.6831345074599439e-25,  1.0000000000000000e+00\n"
"         };\n"
"         ratPolyGate();\n"
"         ForwardEuler(Xr1Gate);\n"
"\n"
"      }\n"
"      //4\n"
"      {\n"
"         const int Mhu_l = 1;\n"
"         const int Mhu_m = 10;\n"
"         const double Mhu_a[] = {\n"
"            2.4919224554816769e-02, -1.0132133908609619e-03,  2.0115370077599340e-05, -2.4787509921554170e-07,  1.9157062546455413e-09, -1.1856625513158431e-11,  4.3373559510489725e-14,  2.9410875655984132e-15, -2.8535346738989432e-17, -2.7319284078381429e-19,  1.0000000000000000e+00\n"
"         };\n"
"         const int Tau_l = 1;\n"
"         const int Tau_m = 10;\n"
"         const double Tau_a[] = {\n"
"            3.2797890377243832e-01,  1.1724387906430438e-06,  3.7188870926079577e-05, -3.3513027526302832e-09,  7.4502277507770199e-09,  1.2888652361968656e-12,  8.0734422444730268e-13,  8.1028416357975461e-16, -1.9648667839939542e-18, -3.5560477657149278e-19,  1.0000000000000000e+00\n"
"         };\n"
"         ratPolyGate();\n"
"         ForwardEuler(Xr2Gate);\n"
"\n"
"      }\n"
"      //5\n"
"      {\n"
"         const int Mhu_l = 5;\n"
"         const int Mhu_m = 5;\n"
"         const double Mhu_a[] = {\n"
"            5.8834144047581471e-01,  1.8554918971499884e-02,  2.3769494810239468e-04,  1.4558400579866515e-06,  3.5608483256093079e-09,  1.0000000000000000e+00,  2.1263634327928394e-03,  5.2685461356465763e-04, -1.2603117542017038e-06,  1.5082411859988069e-08\n"
"         };\n"
"         const int Tau_l = 6;\n"
"         const int Tau_m = 9;\n"
"         const double Tau_a[] = {\n"
"            1.2782667970990670e-03, -5.3480687790940232e-05,  4.1463175377037539e-06, -7.9427267607107334e-08,  2.3181534592617566e-09, -1.8597152793072996e-11,  3.4109130605190490e-13,  2.0323839081334259e-15,  4.3130861032059193e-18,  1.0000000000000000e+00,  4.9865573111495263e-03,  1.0343375562869115e-03,  1.3748473468950345e-05,  3.3217061243887820e-07, -2.7950790321181472e-09\n"
"         };\n"
"         ratPolyGate();\n"
"         ForwardEuler(XsGate);\n"
"      }\n"
"      //6\n"
"      {\n"
"         const int Mhu_l = 5;\n"
"         const int Mhu_m = 8;\n"
"         const double Mhu_a[] = {\n"
"            3.4450952356981063e-02,  3.2569778153723060e-03,  1.4015137277741737e-04,  3.5177311829244476e-06,  5.4744904248961321e-08,  5.1978297648799841e-10,  2.7509401026136440e-12,  6.2017518287576868e-15,  1.0000000000000000e+00, -6.6423744780337163e-02,  2.2756546586311602e-03, -3.0510539067819270e-05,  3.1936830897476799e-07\n"
"         };\n"
"         const int Tau_l = 5;\n"
"         const int Tau_m = 7;\n"
"         const double Tau_a[] = {\n"
"            2.1250755813646321e-01,  5.4375121806696314e-03,  1.2855334968507940e-04,  1.8435064622728614e-06,  2.2201936942212559e-08,  1.6208374761956396e-10,  7.0475957646256506e-13,  1.0000000000000000e+00, -1.1298206006665951e-02,  1.9108836371246739e-05,  3.8822652158898901e-06,  2.4266538879422829e-08\n"
"         };\n"
"         ratPolyGate();\n"
"         ForwardEuler(rGate);\n"
"      }\n"
"      //7\n"
"      {\n"
"         const int Mhu_l = 5;\n"
"         const int Mhu_m = 7;\n"
"         const double Mhu_a[] = {\n"
"            7.4398777641576486e-01,  4.4326834369266523e-02,  1.1359075050132548e-03,  1.5872707499276087e-05,  1.2622142062210935e-07,  5.3559505029241141e-10,  9.3621878808126576e-13,  1.0000000000000000e+00,  2.5428070655441278e-02,  1.7701065218545982e-03,  3.7348639823389864e-06,  2.5522660570009651e-07\n"
"         };\n"
"         const int Tau_l = 7;\n"
"         const int Tau_m = 10;\n"
"         const double Tau_a[] = {\n"
"            1.5068318180921829e+00,  7.8993621885004611e-02,  6.9687982613764387e-03,  1.1019128271535092e-04,  2.3993512027253813e-06,  1.1356515133023819e-08,  3.1836750643625883e-10, -1.6718952430054101e-12,  2.5579267994953318e-14,  1.0819050682737564e-16,  1.0000000000000000e+00, -6.8101965713816826e-02,  3.7430651166330204e-03,  6.6379282851072469e-07,  9.3098141480209233e-07,  3.2614589730112466e-08,  4.4763500756297451e-10\n"
"         };\n"
"         ratPolyGate();\n"
"         ForwardEuler(dGate);\n"
"      }\n"
"      //8\n"
"      {\n"
"         const int Mhu_l = 5;\n"
"         const int Mhu_m = 8;\n"
"         const double Mhu_a[] = {\n"
"            5.4304780098652715e-02, -4.4962960006268609e-03,  1.7191393153273102e-04, -3.9003819270270578e-06,  5.4530196576796994e-08, -4.1503002319748868e-10,  8.3481595187685469e-13,  5.8904944849510583e-15,  1.0000000000000000e+00,  5.2318145763546209e-02,  1.6278824780069224e-03,  1.6010969045398731e-05,  1.9936720027259198e-07\n"
"         };\n"
"         const int Tau_l = 13;\n"
"         const int Tau_m = 11;\n"
"         const double Tau_a[] = {\n"
"            8.7301971950165709e-03,  1.0619909606051185e-03,  7.4439522770957243e-05,  3.5189522703546544e-06,  1.2763275154107722e-07,  3.9329666177349803e-09,  1.0963505894631244e-10,  2.5405559223755679e-12,  4.1926387533661996e-14,  4.1426318225610267e-16,  1.8370569627175333e-18,  1.0000000000000000e+00,  5.3496360858059103e-02,  1.0576677271235941e-02,  6.0292981147301659e-04,  2.4470545201049391e-05,  8.8342556285406000e-07,  2.5756751722220191e-08,  5.6840520101297543e-10,  8.9563864907355775e-12,  8.8109045583454404e-14,  4.1790988048136637e-16,  3.0481930918230229e-19,  7.4195798494855779e-22\n"
"         };\n"
"         ratPolyGate();\n"
"         ForwardEuler(fGate);\n"
"      }\n"
"      //9\n"
"      {\n"
"         const int Mhu_l = 8;\n"
"         const int Mhu_m = 5;\n"
"         const double Mhu_a[] = {\n"
"            3.3448475117101473e-01,  1.7908684328425843e-02,  4.6473146231373311e-04,  4.9947394623057196e-06,  4.6013519663601120e-08,  1.0000000000000000e+00,  5.5442666103601546e-02,  1.3608399868819414e-03,  1.6306821023781025e-05,  1.2063672477786111e-07,  1.7672562299303354e-10, -5.4528541905004782e-13, -2.8928059365721159e-15\n"
"         };\n"
"         const int Tau_l = 11;\n"
"         const int Tau_m = 12;\n"
"         const double Tau_a[] = {\n"
"            3.0215260344315119e-02,  2.9905565403315244e-03,  1.8698675038604531e-04,  7.8302945745729333e-06,  2.4353851952274140e-07,  5.9716586064763031e-09,  1.4320157611618665e-10,  3.5798223309278773e-12,  7.3139244140308910e-14,  9.1797821675112518e-16,  5.5127688299509001e-18,  6.0160642283988931e-21,  1.0000000000000000e+00, -8.8631007579573415e-02,  5.6222208830099462e-03,  9.0331117197141856e-05,  8.5481830250918979e-07,  2.2780719183249479e-07,  7.0378379382184305e-09,  1.0211239994063565e-10,  1.8042686984162349e-12,  2.9931381582699879e-14,  2.1128053668752277e-16\n"
"         };\n"
"         ratPolyGate();\n"
"         ForwardEuler(f2Gate);\n"
"      }\n"
"      //10\n"
"      {\n"
"         const int Mhu_l = 6;\n"
"         const int Mhu_m = 12;\n"
"         const double Mhu_a[] = {\n"
"            3.6932725836170428e-09,  2.5736169972850485e-10, -7.5177270851662063e-11, -1.7783383158540549e-12,  2.3780959584764758e-13,  3.9995914648833690e-15, -2.6002507057823451e-16, -4.0173750386716292e-18,  1.1567747025479995e-19,  1.8190095403795220e-21, -1.7772014587013272e-23, -3.0513026345262148e-25,  1.0000000000000000e+00,  4.7761481528554264e-02,  9.1744037118706855e-04,  8.8560239615590481e-06,  4.2943755440989806e-08,  8.3656852902593096e-11\n"
"         };\n"
"         const int Tau_l = 1;\n"
"         const int Tau_m = 1;\n"
"         const double Tau_a[] = {\n"
"            1.4925373134328358e-03,  1.0000000000000000e+00\n"
"         };\n"
"         ratPolyGate();\n"
"         ForwardEuler(jLGate);\n"
"      }\n"
"      //11\n"
"      {\n"
"         const int Mhu_l = 7;\n"
"         const int Mhu_m = 8;\n"
"         const double Mhu_a[] = {\n"
"            3.6866939888900814e-03, -4.5574539113452475e-04,  2.6788146214637530e-05, -9.9102133439606097e-07,  2.5613910193807883e-08, -4.6400777680998212e-10,  5.2352368345974410e-12, -2.6731132916017718e-14,  1.0000000000000000e+00,  7.5601232522433295e-02,  2.5744765696710309e-03,  4.7389711229137498e-05,  5.8624244628409614e-07,  3.4873125981144344e-09,  2.0871941515527373e-11\n"
"         };\n"
"         const int Tau_l = 6;\n"
"         const int Tau_m = 7;\n"
"         const double Tau_a[] = {\n"
"            5.2008879668973697e-02,  2.9402318950257256e-03,  9.6193142778593343e-05,  1.9694354926101320e-06,  2.5766968301605511e-08,  1.9309877695248136e-10,  6.5471492087484267e-13,  1.0000000000000000e+00, -2.1715782138805460e-02,  2.0915626226951107e-03,  1.8073629807905746e-07,  2.1067307291025196e-07,  2.7434857789930993e-09\n"
"         };\n"
"         ratPolyGate();\n"
"         ForwardEuler(sGate);\n"
"      }   \n"
"}\n"
"\n"
;

void actualCalcJitify(const double dt, const int nCells_, const double Vm[], const double iStim[], double dVm[], double* stateData) {
   const int blockSize=512;

/*
   static jitify::JitCache kernel_cache;
   jitify::Program program = kernel_cache.program(SimpleTT06_ker,0);

   program
      .kernel("actualCalcCudaJitify")
      .instantiate()
      .configure(dim3((nCells_+blockSize-1)/blockSize,1,1),dim3(blockSize,1,1))
      .launch(dt,nCells_,Vm,iStim,dVm,stateData);
*/
}


   
}
