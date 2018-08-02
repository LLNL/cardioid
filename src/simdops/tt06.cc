

#include <new>
#include <limits>
#include <cstdlib>
#include <cstddef>
#include <cmath>
#include <vector>
#include <iostream>
#include <sys/time.h>
#include <cstring>
#include <stdint.h>
#include "AlignedAllocator.hh"
#include "simdops/simdops.hpp"

extern "C" {
void HPM_Init(); 
void HPM_Start(char * this_label);
void HPM_Stop(char * this_label);
void  HPM_Print();
}

#if 0
#include <omp.h>
#else
int omp_get_max_threads() {return 1; }
int omp_get_num_procs() { return 1; }
int omp_get_thread_num() { return 0; }
#endif

using namespace std;
using namespace simdops;

typedef double Real;
typedef simdops::float64v VReal;

typedef VReal VRealArray;
#define mfarraydef(x) make_float(x)
#define mfarrayuse(x) (x)
#define RESOLVE(x) (x)
inline Real toDouble(const VReal x) {
   Real SIMDOPS_ALIGN(SIMDOPS_FLOAT64V_WIDTH*sizeof(double)) el[SIMDOPS_FLOAT64V_WIDTH];
   store(el, x);
   return el[0];
}

#define vectorSimd(T) vector<T, AlignedAllocator<T, SIMDOPS_FLOAT64V_WIDTH*sizeof(double)> >



inline uint64_t getTime() {
  return 0;
}

struct Timer
{
   uint64_t start_;

   inline void tic() {
      start_ = getTime();
   }
   inline uint64_t toc() {
      return getTime()-start_;
   }
};

struct State 
{
   //EDIT_STATE
   VReal *f2Gate;
   VReal *fGate;
   VReal *dGate;
   VReal *mGate;
   VReal *jGate;
   VReal *hGate;
   VReal *rGate;
   VReal *sGate;
   VReal *Xr1Gate;
   VReal *Xr2Gate;
   VReal *XsGate;
   VReal *jLGate;
   VReal *Na_i;
   VReal *Ca_i;
   VReal *Ca_ss;
   VReal *Ca_sr;
   VReal *fCass;
   VReal *dVK_i;
   VReal *R_prime;
} SIMDOPS_ALIGN(SIMDOPS_FLOAT64V_WIDTH*sizeof(double));

#define CUBE(x) ((x)*(x)*(x))
#define SQ(x) ((x)*(x))
#define sigm(x)   div((x),1+(x))
#define logSeries(x)    (log(1+(x)) )

typedef Real* alignedRealPtr ;//__attribute__((align_value(SIMD_WIDTH)));
typedef State* alignedStatePtr ;//__attribute__((align_value(SIMD_WIDTH)));
typedef VReal* alignedVRealPtr; //__attribute__((align_value(SIMD_WIDTH)));

uint64_t globalNonGate;
uint64_t globalGate;

void nonGates(const Real _dt, const int nSimd, const int simdBegin, const alignedRealPtr Vm, const alignedRealPtr iStim, alignedRealPtr dVm, alignedStatePtr state_)
{

  const VReal dt = splat(&_dt);
   for (int ii=simdBegin; ii<nSimd+simdBegin; ++ii)
   {
      //set per-cell flags
      //EDIT_PERCELL_FLAGS
      
      //set per-cell parameters
      //EDIT_PERCELL_PARAMETERS

      const Real pcnst_0 = 8314.472;
      const Real pcnst_1 = 310;
      const Real pcnst_2 = 96485.3415;
      const Real pcnst_9 = 0.03;
      const Real pcnst_10 = 5.4;
      const Real pcnst_11 = 140;
      const Real pcnst_12 = 2;
      const Real _c1 = pcnst_2/(pcnst_0*pcnst_1);
      const VReal c1 = make_float(_c1 );
      const Real _c2 = pcnst_9;
      const VReal c2 = make_float(_c2 );
      const Real _c3 = -1/_c1;
      const VReal c3 = make_float(_c3 );

      //const Real _c4 = -c3*log(pcnst_11);
      //const Real _c5 = -c3*log(pcnst_10);
      //const Real _c6 = -c3*log(pcnst_10+pcnst_9*pcnst_11);
      //const Real _c8 = -0.5*c3*log(pcnst_12);
      //if (0) {
      //  printf("%.17g %.17g %.17g %.17g\n", c4,c5,c6,c8);
      //}

      const Real _c4 = 132.00985294338352;
      const VReal c4 = make_float(_c4 );
      const Real _c5 = 45.050058022436289;
      const VReal c5 = make_float(_c5 );
      const Real _c6 = 60.420198081560486;
      const VReal c6 = make_float(_c6 );
      const Real _c8 = 9.2582839417106104;
      const VReal c8 = make_float(_c8 );
      const Real _P_NaK=2.724;
      const VReal P_NaK= make_float(_P_NaK);
      const Real _g_Ks = 0.392;
      const VReal g_Ks = make_float(_g_Ks );
      const Real _g_Kr = 0.153;
      const VReal g_Kr = make_float(_g_Kr );
      const Real _g_to = 0.073;
      const VReal g_to = make_float(_g_to );
      const Real _g_NaL=0.0;
      const VReal g_NaL= make_float(_g_NaL);

      const Real pcnst_3 = 0.185;
      const Real pcnst_4 = 0.016404;
      const Real pcnst_5 = 10;
      const Real pcnst_6 = 1000;
      const Real pcnst_7 = 1;
      const Real pcnst_8 = 52;
      const Real pcnst_13 = 5.405;
      const Real pcnst_14 = 0.153;
      const Real pcnst_16 = 14.838;
      const Real pcnst_17 = 0.00029;
      const Real pcnst_18 = 0.0000398;
      const Real pcnst_19 = 0.000592;
      const Real pcnst_21 = 2.724;
      const Real pcnst_22 = 1;
      const Real pcnst_23 = 40;
      const Real pcnst_24 = 1000;
      const Real pcnst_25 = 0.1;
      const Real pcnst_26 = 2.5;
      const Real pcnst_27 = 0.35;
      const Real pcnst_28 = 1.38;
      const Real pcnst_29 = 87.5;
      const Real pcnst_30 = 0.1238;
      const Real pcnst_31 = 0.0005;
      const Real pcnst_32 = 0.0146;
      const Real pcnst_33 = 0.15;
      const Real pcnst_34 = 0.045;
      const Real pcnst_35 = 0.06;
      const Real pcnst_36 = 0.005;
      const Real pcnst_37 = 1.5;
      const Real pcnst_38 = 2.5;
      const Real pcnst_39 = 1;
      const Real pcnst_40 = 0.102;
      const Real pcnst_41 = 0.0038;
      const Real pcnst_42 = 0.00025;
      const Real pcnst_43 = 0.00036;
      const Real pcnst_44 = 0.006375;
      const Real pcnst_45 = 0.2;
      const Real pcnst_46 = 0.001;
      const Real pcnst_47 = 10;
      const Real pcnst_48 = 0.3;
      const Real pcnst_49 = 0.4;
      const Real pcnst_50 = 0.00025;
      const Real pcnst_51 = 0.001094;
      const Real pcnst_52 = 0.00005468;

      const Real _c9 = -pcnst_3/(pcnst_4*pcnst_2);
      const VReal c9 = make_float(_c9 );


      const Real _c7 =  pcnst_19;
      const VReal c7 = make_float(_c7 );
      //c11= pcsnt[14]*sqrt(pcnst_10/5.4);
      //const Real _c11= sqrt(pcnst_10/5.4);
      const Real _c11 = 1;
      const VReal c11 = make_float(_c11 );
      //const Real _c12= pcnst_13*sqrt(pcnst_10/5.4);
      const Real _c12 = 5.4050000000000002;
      const VReal c12 = make_float(_c12 );
      const Real _c13= pcnst_3/(2.0*pcnst_52*pcnst_2*_c9);
      const VReal c13= make_float(_c13);
      const Real _c14 = pcnst_51/pcnst_52;
      const VReal c14 = make_float(_c14 );
      const Real _c15 = -pcnst_52/pcnst_4;
      const VReal c15 = make_float(_c15 );
      const Real _c16 = pcnst_51/pcnst_4;
      const VReal c16 = make_float(_c16 );
      const Real _c17 = pcnst_35/(pcnst_33*pcnst_34)*_c9;
      const VReal c17 = make_float(_c17 );
      const Real _c18 = pcnst_34*pcnst_38/_c9 ;
      const VReal c18 = make_float(_c18 );
      const Real _c19  = -pcnst_34*(pcnst_38-pcnst_39)/_c9;
      const VReal c19  = make_float(_c19  );
      const Real _c20  = pcnst_16;
      const VReal c20  = make_float(_c20  );
      const Real _c21  = pcnst_17;
      const VReal c21  = make_float(_c21  );
      const Real _c22  = 1/_c9;
      const VReal c22  = make_float(_c22  );
      const Real _c23  = pcnst_41/(_c15*_c9);
      const VReal c23  = make_float(_c23  );
      const Real _c24  =  pcnst_30;
      const VReal c24  = make_float(_c24  );
      const Real _c25  =  1.0/pcnst_23;
      const VReal c25  = make_float(_c25  );
      const Real _c26  =  1.0/pcnst_31;
      const VReal c26  = make_float(_c26  );
      const Real _c27  =  1.0/pcnst_42;
      const VReal c27  = make_float(_c27  );
      //const Real c28  =  1.0/sqrt(pcnst_45*pcnst_46);
      const Real _c28 = 70.710678118654755;
      const VReal c28 = make_float(_c28 );
      const Real _c29  =  pcnst_46*_c28;
      const VReal c29  = make_float(_c29  );
      const Real _c30  =  1.0/pcnst_37;
      const VReal c30  = make_float(_c30  );
      //const Real _c31  =  1.0/sqrt(pcnst_47*pcnst_48);
      const Real _c31 = 0.57735026918962584;
      const VReal c31 = make_float(_c31 );
      const Real _c32  =  pcnst_48*_c31;
      const VReal c32  = make_float(_c32  );
      //const Real _c33  =  1.0/sqrt(pcnst_49*pcnst_50);
      const Real _c33 = 100;
      const VReal c33 = make_float(_c33 );
      const Real _c34  =  pcnst_50*_c33;
      const VReal c34  = make_float(_c34  );
      const Real _c36  =  pcnst_36/_c9;
      const VReal c36  = make_float(_c36  );
      const Real _c40  =  pcnst_40/_c9;
      const VReal c40  = make_float(_c40  );
      const Real _c43  =  pcnst_43/_c9;
      const VReal c43  = make_float(_c43  );
      const Real _c44  =  pcnst_44/_c9;
      const VReal c44  = make_float(_c44  );

      const Real _f1 = _c1;
      const VReal f1 = make_float(_f1 );
      const Real _f2 =  -2.0*pcnst_10/(pcnst_10+pcnst_22);
      const VReal f2 = make_float(_f2 );
      const Real _f3 =  ((CUBE(pcnst_29)+CUBE(pcnst_11))*(pcnst_28+pcnst_12))/(pcnst_24*pcnst_12);
      const VReal f3 = make_float(_f3 );
      const Real _f4 =  _f3*pcnst_25;
      const VReal f4 = make_float(_f4 );
   
      const Real _f5 =  pcnst_27*_f1;
      const VReal f5 = make_float(_f5 );
      const Real _f6 =  (CUBE(pcnst_11)*pcnst_26/pcnst_12);
      const VReal f6 = make_float(_f6 );
      const Real _f7 = pcnst_18*pcnst_2*_f1;
      const VReal f7 = make_float(_f7 );
      const Real _f7a = 0.5*pcnst_18*pcnst_2;
      const VReal f7a = make_float(_f7a );
      const Real _f9  = 4.0*pcnst_12;
      const VReal f9  = make_float(_f9  );
      const Real _f9a = 4.0*pcnst_12*_f7a;
      const VReal f9a = make_float(_f9a );
      const Real _f10 = pcnst_32;
      const VReal f10 = make_float(_f10 );

      const Real  _c_K1   =  1;
      const VReal c_K1   = make_float(_c_K1   );
      const Real  _c_Na   =  1 *_c20;
      const VReal c_Na   = make_float(_c_Na   );
      const Real  _c_bNa  =  1*_c21;
      const VReal c_bNa  = make_float(_c_bNa  );
      const Real  _c_CaL  =  1;
      const VReal c_CaL  = make_float(_c_CaL  );
      const Real  _c_bCa  =  1*_c7;
      const VReal c_bCa  = make_float(_c_bCa  );
      const Real  _c_NaCa =  1;
      const VReal c_NaCa = make_float(_c_NaCa );
      const Real  _c_pCa  =  1*_c24;
      const VReal c_pCa  = make_float(_c_pCa  );
      const Real  _c_pK   =  1;
      const VReal c_pK   = make_float(_c_pK   );
      const Real  _c_up   = 1*_c44;
      const VReal c_up   = make_float(_c_up   );
      const Real  _c_leak = 1*_c43;
      const VReal c_leak = make_float(_c_leak );
      const Real  _c_xfer = 1*_c23;
      const VReal c_xfer = make_float(_c_xfer );
      const Real  _c_rel  = 1*_c40;
      const VReal c_rel  = make_float(_c_rel  );

      const Real _c_NaK = _P_NaK;
      const VReal c_NaK = make_float(_c_NaK );
      const Real _c_Ks  = _g_Ks;
      const VReal c_Ks  = make_float(_c_Ks  );
      const Real _c_Kr  = _g_Kr;
      const VReal c_Kr  = make_float(_c_Kr  );
      const Real _c_to  = _g_to;
      const VReal c_to  = make_float(_c_to  );
      const Real _c_NaL = _g_NaL;
      const VReal c_NaL = make_float(_c_NaL );

      //set Vm
      const VReal thisVm = load(&Vm[ii*SIMDOPS_FLOAT64V_WIDTH]);
      const VReal istim = load(&iStim[ii*SIMDOPS_FLOAT64V_WIDTH]);

      //set all state variables
      //EDIT_STATE
      const VReal f2Gate=state_->f2Gate[ii];
      const VReal fGate=state_->fGate[ii];
      const VReal dGate=state_->dGate[ii];
      const VReal mGate=state_->mGate[ii];
      const VReal jGate=state_->jGate[ii];
      const VReal hGate=state_->hGate[ii];
      const VReal rGate=state_->rGate[ii];
      const VReal sGate=state_->sGate[ii];
      const VReal Xr1Gate=state_->Xr1Gate[ii];
      const VReal Xr2Gate=state_->Xr2Gate[ii];
      const VReal XsGate=state_->XsGate[ii];
      const VReal jLGate=state_->jLGate[ii];
      const VReal _Na_i=state_->Na_i[ii];
      const VReal _Ca_i=state_->Ca_i[ii];
      const VReal _Ca_ss=state_->Ca_ss[ii];
      const VReal _Ca_SR=state_->Ca_sr[ii];
      const VReal _fCass=state_->fCass[ii];
      const VReal _dVK_i=state_->dVK_i[ii];
      const VReal _R_prime=state_->R_prime[ii];

      //Real *states = cell[ii].state;
      VReal dVR = make_float(0.0);
      //Real itmp0,itmp5,itmp6;
      VReal I_K1,I_Kr,I_Ks,I_Na,I_bNa,I_CaL,I_bCa,I_to,I_NaK,I_NaCa,I_pCa,I_pK,I_NaL,I_leak,I_up,I_rel,I_xfer;
      VReal I_sum,I_delta;


      //  Update Ca concentration;
      {
         VReal fv1; {
            const VReal x = thisVm;
            const int l=1,m=6;
            const VRealArray a[]={mfarraydef(1.5757133041269722e-04),mfarraydef(2.4132516902459823e-06),mfarraydef(1.4605316459098998e-08),mfarraydef(6.0427222047915772e-11),mfarraydef(2.6019285647913959e-13),mfarraydef(6.0217000704020008e-16),mfarraydef(1.0000000000000000e+00)
            };
            VReal sum1 = make_float(0);
            #pragma unroll(m)
            for (int j = m-1; j >= 0; j--)
              sum1 = add(mfarrayuse(a[j]),mul(x,sum1));

            fv1 = sum1;
         }

         VReal fv2; {
            const VReal x = thisVm;
            const int l=1,m=8;
            const VRealArray a[]={mfarraydef(5.4050480762039217e+02),mfarraydef(-1.1953189161096423e+01),mfarraydef(1.1865349408070526e-01),mfarraydef(-6.0092610681801851e-04),mfarraydef(4.0077214280113720e-07),mfarraydef(2.9421099617665475e-08),mfarraydef(-2.5528046472894199e-10),mfarraydef(-1.5237049656618966e-12),mfarraydef(1.0000000000000000e+00)
            };
            VReal sum1 = make_float(0);
            #pragma unroll(m)
            for (int j = m-1; j >= 0; j--)
               sum1 = add(mfarrayuse(a[j]),mul(x,sum1));

            fv2 = sum1;
         }
         VReal x1 = _Ca_i * c26;
         VReal x2 = SQ(_Ca_i * c27);
         VReal x3 = SQ(_Ca_i * c28 + c29);
         VReal sigm1 = sigm(x1);
         VReal sigm2 = sigm(x2);
         VReal sigm3 = sigm(x3);
         VReal dV3 = thisVm - 0.5 * c3 * log(_Ca_i) - c8;

         I_NaCa = c_NaCa * (CUBE(_Na_i) * fv1 - _Ca_i * fv2);

         I_bCa = c_bCa * dV3;
         I_pCa = c_pCa * sigm1;
         I_up  =  c44 * sigm2;
         I_leak=  c43 * (_Ca_SR - _Ca_i);
         I_xfer = c23 * (_Ca_ss - _Ca_i);
         I_delta = I_leak - I_up; // I_detal = -itmp5
         I_sum =   I_bCa + I_pCa; // itmp4 =I_bCa+I_pCa;
         state_->Ca_i[ii]   = _Ca_i + (dt * c9) * (sigm3 * (0.5 * I_sum - I_NaCa + I_xfer * c15 + I_delta*c16));
         dVR  = dVR - I_sum;
         //if (ii %4 ==0) printf("\n%d dVR=%14.12f ", ii,dVR);
         //else printf("%14.12f ", ii,I_sum);
      }

      VReal iK;
      //  Update K and Na concentrations;
      {
         VReal _K_i  = c9 * (_dVK_i - thisVm);
         VReal x0 = _Na_i * c25;

         VReal dV0 = thisVm - c3*log(_K_i) - c5;
         VReal dV1 = thisVm - c3*log(_Na_i) - c4;
         VReal dV2 = dV0 - c3 * logSeries(RESOLVE(c2 * _Na_i/_K_i)) + c5 - c6; // Assumption:  c2*_Na_i/_K_i is small;

         VReal fv0; {
            const VReal x = thisVm;
            const int l=1,m=10;
            const VRealArray a[]={mfarraydef(-1.4549948900403902e+00),mfarraydef(-2.2427240624934000e-03),mfarraydef(2.8701981811073087e-05),mfarraydef(-2.9352599282738483e-07),mfarraydef(1.8687813196724116e-09),mfarraydef(-6.1533240726750350e-12),mfarraydef(-3.0715514390797755e-14),mfarraydef(2.1068763323371155e-15),mfarraydef(-1.5741479441621198e-17),mfarraydef(-1.6908987107052300e-19),mfarraydef( 1.0000000000000000e+00)
            };
            VReal sum1 = make_float(0);
            #pragma unroll(m)
            for (int j = m-1; j >= 0; j--)
               sum1 = add(mfarrayuse(a[j]),mul(x,sum1));

            fv0 = sum1;

         }
         VReal fv5; {
            const VReal x = thisVm;
            const int l=5,m=8;
            const VRealArray a[]={mfarraydef(2.1983170611730303e-04),mfarraydef(2.1698309464995643e-05),mfarraydef(9.8317099143409939e-07),mfarraydef(2.6134710139327942e-08),mfarraydef(4.3109819211796886e-10),mfarraydef(4.3229574960011243e-12),mfarraydef(2.4019472580987271e-14),mfarraydef(5.6472223920064411e-17),mfarraydef(1.0000000000000000e+00),mfarraydef(-6.6157324225412789e-02),mfarraydef(2.0076473883411361e-03),mfarraydef(-2.6510671504853674e-05),mfarraydef(2.2430697618084790e-07)
            };
            VReal sum1 = make_float(0);
            #pragma unroll(m)
            for (int j = m-1; j >= 0; j--)
               sum1 = add(mfarrayuse(a[j]),mul(x,sum1));

            VReal sum2 = make_float(0);
            int k = m + l - 1;
            #pragma unroll(l) 
            for (int j = k; j >= m; j--)
               sum2 = add(mfarrayuse(a[j]),mul(x,sum2));
            fv5 = div(sum1,sum2);
         }

         VReal fv6; {
            const VReal x = dV0;
            const int l=11,m=7;
            const VRealArray a[]={mfarraydef(2.9794200055170916e-01),mfarraydef(-8.8983788808079303e-03),mfarraydef(7.9933972784803960e-03),mfarraydef(-1.7599632267422598e-04),mfarraydef(1.5074733144798214e-06),mfarraydef(-5.9120168317614184e-09),mfarraydef(8.9356468869962987e-12),mfarraydef(1.0000000000000000e+00),mfarraydef(2.1722631714279653e-01),mfarraydef(2.1762877423819233e-02),mfarraydef(1.5923509782789353e-03),mfarraydef(-7.4897050385073499e-05),mfarraydef(3.3873719083779318e-06),mfarraydef(-5.3689489360972982e-08),mfarraydef(8.8474419958442444e-10),mfarraydef(-4.4602715802630215e-12),mfarraydef(2.0785335146560804e-14),mfarraydef(1.9580847130746169e-16)
            };
            VReal sum1 = make_float(0);
            #pragma unroll(m)
            for (int j = m-1; j >= 0; j--)
               sum1 = add(mfarrayuse(a[j]),mul(x,sum1));

            VReal sum2 = make_float(0);
            int k = m + l - 1;
            #pragma unroll(l) 
            for (int j = k; j >= m; j--)
               sum2 = add(mfarrayuse(a[j]),mul(x,sum2));
            fv6 = div(sum1,sum2);
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
         
         VReal iNa =  3 * I_NaCa - 1.5 * I_NaK + I_Na + I_bNa + I_NaL;
         iK =  I_Ks + I_NaK + I_pK + I_K1 + I_to + I_Kr;


         state_->dVK_i[ii] = state_->dVK_i[ii] + dt * iK;
         dVR    =  dVR - (iNa + iK - 2.0 * I_NaCa);
         state_->Na_i[ii] = _Na_i + (dt * c9) * iNa;
      }

      //  Update Ca_SS, Ca_SR, R_prime concentrations and fCass gate;
      {
         VReal fv3; {
            const VReal x = thisVm;
            const int l=1,m=11;
            const VRealArray a[]={mfarraydef(1.0396540551771880e+00),mfarraydef(4.6058010729117471e-02),mfarraydef(7.9107591823998407e-04),mfarraydef(4.3024156297446584e-06),mfarraydef(-4.9702295945293119e-08),mfarraydef(-6.9663904851401655e-10),mfarraydef(2.4991578397134279e-12),mfarraydef(7.8166616197552177e-14),mfarraydef(1.3685444078468510e-16),mfarraydef(-3.6584729114311930e-18),mfarraydef(-1.7450648728130940e-20),mfarraydef(1.0000000000000000e+00)
            };
            VReal sum1 = make_float(0);
            #pragma unroll(m)
            for (int j = m-1; j >= 0; j--)
               sum1 = add(mfarrayuse(a[j]),mul(x,sum1));

            fv3 = sum1;
         }

         VReal fv4; {
            const VReal x = thisVm;
            const int l=1,m=9;
            const VRealArray a[]={mfarraydef(2.5568410537089715e+01),mfarraydef(-7.8144592191919016e-01),mfarraydef(6.3170153846905984e-03),mfarraydef(3.4110729724745464e-05),mfarraydef(-3.7991861595738163e-07),mfarraydef(-5.2236184692414528e-09),mfarraydef(1.1857805262424718e-11),mfarraydef(4.5506699371915196e-13),mfarraydef(1.8913634488808377e-15),mfarraydef(1.0000000000000000e+00)
            };
            VReal sum1 = make_float(0);
            #pragma unroll(m)
            for (int j = m-1; j >= 0; j--)
               sum1 = add(mfarrayuse(a[j]),mul(x,sum1));

            fv4 = sum1;
         }
         VReal x4 = SQ(_Ca_SR * c30);
         VReal x5 = SQ(_Ca_SR * c31 + c32);
         VReal x6 = SQ(_Ca_ss * c33 + c34);
         VReal sigm4 = sigm(x4);
         VReal sigm5 = sigm(x5);
         VReal sigm6 = sigm(x6);
         VReal tmp8  = (c18 + c19 * sigm4); //Sigm4
         VReal tmp9  = tmp8 * _Ca_ss + c36;

         I_CaL = c_CaL * dGate * fGate * f2Gate * _fCass * (fv4 - _Ca_ss * fv3);  // renamed itmp1 to I_CaL

         VReal O = (SQ(_Ca_ss) * _R_prime)/((tmp8 * c17 + SQ(_Ca_ss)));
         I_rel =c40 * O * (_Ca_SR - _Ca_ss);

         state_->Ca_ss[ii]   = _Ca_ss   + (dt * c9) * sigm6 * (I_xfer + I_rel * c14 + I_CaL * c13);
         state_->Ca_sr[ii]   = _Ca_SR   - (dt * c9) * sigm5 * (I_delta + I_rel);
         state_->R_prime[ii] = _R_prime + (dt * c9) * (c36 - tmp9 * _R_prime);

         //#if fCassForm == TT06
         VReal t1 = (1.0)/((1.0 + SQ(20 * _Ca_ss)));
         VReal mhu = 0.600000 * t1 + 0.4000000;
         VReal tauR = (1.0)/((80.0 * t1 + 2.0));
         //#endif

         /*#if  fCassForm == RICE
         Real xCa_ss = 1000*_Ca_ss;
         Real xCa_i  = 1000*_Ca_i ;
         Real mhu    = 0.6/(1.0 + xCa_ss * sqrt(xCa_ss) + xCa_i * sqrt(xCa_i)) + 0.4;
         Real tauR   = 0.005/mhu;
         #endif
         */
         state_->fCass[ii]   = _fCass   + dt*(mhu - _fCass) * tauR;

         dVR = dVR + I_CaL;
      }

      store(&dVm[ii*SIMDOPS_FLOAT64V_WIDTH], dVR);

   }
}

#define RushLarsen(name) thisState_[ii] = thisState_[ii] - (mhu - thisState_[ii])*expm1(RESOLVE(-1*tauR));
#define ForwardEuler(name) thisState_[ii] = thisState_[ii] + (mhu - thisState_[ii])*tauR

void mGate(const Real _dt, const int nSimd, const int simdBegin, const alignedRealPtr Vm, const alignedRealPtr iStim, alignedRealPtr dVm, alignedVRealPtr thisState_)
{

  const VReal dt = splat(&_dt);
   for (int ii=simdBegin; ii<nSimd+simdBegin; ++ii)
   {
      const VReal thisVm = load(&Vm[ii*SIMDOPS_FLOAT64V_WIDTH]);

      //0
      {
      const int Mhu_l = 10;
      const int Mhu_m = 5;
      const VRealArray Mhu_a[] = { mfarraydef(9.9632117206253790e-01),mfarraydef(  4.0825738726469545e-02),mfarraydef(  6.3401613233199589e-04),mfarraydef(  4.4158436861700431e-06),mfarraydef(  1.1622058324043520e-08),mfarraydef(  1.0000000000000000e+00),mfarraydef(  4.0568375699663400e-02),mfarraydef(  6.4216825832642788e-04),mfarraydef(  4.2661664422410096e-06),mfarraydef(  1.3559930396321903e-08),mfarraydef( -1.3573468728873069e-11),mfarraydef( -4.2594802366702580e-13),mfarraydef(  7.6779952208246166e-15),mfarraydef(  1.4260675804433780e-16),mfarraydef( -2.6656212072499249e-18)};
      const int Tau_l = 1;
      const int Tau_m = 18;
      const VRealArray Tau_a[] = {mfarraydef(1.7765862602413648e+01),mfarraydef(  5.0010202770602419e-02),mfarraydef( -7.8002064070783474e-04),mfarraydef( -6.9399661775931530e-05),mfarraydef(  1.6936588308244311e-06),mfarraydef(  5.4629017090963798e-07),mfarraydef( -1.3805420990037933e-08),mfarraydef( -8.0678945216155694e-10),mfarraydef(  1.6209833004622630e-11),mfarraydef(  6.5130101230170358e-13),mfarraydef( -6.9931705949674988e-15),mfarraydef( -3.1161210504114690e-16),mfarraydef(  5.0166191902609083e-19),mfarraydef(  7.8608831661430381e-20),mfarraydef(  4.3936315597226053e-22),mfarraydef( -7.0535966258003289e-24),mfarraydef( -9.0473475495087118e-26),mfarraydef( -2.9878427692323621e-28),mfarraydef(  1.0000000000000000e+00)};

      VReal sum1,sum2;
      const VReal x = thisVm;

      sum1 = make_float(0);
      #pragma unroll(Mhu_m)
      for (int j = Mhu_m-1; j >= 0; j--)
         sum1 = add(mfarrayuse(Mhu_a[j]),mul(x,sum1));

      sum2 = make_float(0);
      int k = Mhu_m + Mhu_l - 1;
      #pragma unroll(Mhu_l)
      for (int j = k; j >= Mhu_m; j--)
         sum2 = add(mfarrayuse(Mhu_a[j]),mul(x,sum2));
      VReal mhu = div(sum1,sum2);

      sum1 = make_float(0);
      #pragma unroll(Tau_m)
      for (int j = Tau_m-1; j >= 0; j--)
         sum1 = add(mfarrayuse(Tau_a[j]),mul(x,sum1));

      sum2 = make_float(0);
      k = Tau_m + Tau_l - 1;
      #pragma unroll(Tau_l)
      for (int j = k; j >= Tau_m; j--)
         sum2 = add(mfarrayuse(Tau_a[j]),mul(x,sum2));
      VReal tauR = div(sum1,sum2)*dt;
      RushLarsen(mGate);
      }
   }
}
void hGate(const Real _dt, const int nSimd, const int simdBegin, const alignedRealPtr Vm, const alignedRealPtr iStim, alignedRealPtr dVm, alignedVRealPtr thisState_)
{

  const VReal dt = splat(&_dt);
   for (int ii=simdBegin; ii<nSimd+simdBegin; ++ii)
   {
      const VReal thisVm = load(&Vm[ii*SIMDOPS_FLOAT64V_WIDTH]);
      //1
      {
         const int Mhu_l = 7;
         const int Mhu_m = 9;
         const VRealArray Mhu_a[] = {
           mfarraydef(3.1263044718396443e-08),mfarraydef(6.1534783640833612e-09),mfarraydef(-3.7646458226792461e-10),mfarraydef(-4.2451512350227289e-11),mfarraydef(  1.5216916689168596e-12),mfarraydef(  3.9843764909619204e-14),mfarraydef( -7.7509942611798927e-16),mfarraydef( -3.3646392579644030e-17),mfarraydef(  5.7566741989713017e-19),mfarraydef(  1.0000000000000000e+00),mfarraydef(  6.4169311162380532e-02),mfarraydef(  1.7346341071731517e-03),mfarraydef(  2.5211977280197058e-05),mfarraydef(  2.0770353677292379e-07),mfarraydef(  9.2161835229419180e-10),mfarraydef(  1.7368930279808849e-12)
         };
         const int Tau_l = 11;
         const int Tau_m = 11;
         const VRealArray Tau_a[] = {
           mfarraydef(4.2835355918249460e+00),mfarraydef(  6.4588153679793869e-01),mfarraydef(  4.4298327200183743e-02),mfarraydef(  1.8236040756549200e-03),mfarraydef(  5.0201035822685939e-05),mfarraydef(  9.7279248118950919e-07),mfarraydef(  1.3522528070830968e-08),mfarraydef(  1.3356095236104801e-10),mfarraydef(  8.9640565060229970e-13),mfarraydef(  3.6781270948185535e-15),mfarraydef(  6.9724701901279024e-18),mfarraydef(  1.0000000000000000e+00),mfarraydef(  1.2584512299220066e-01),mfarraydef(  7.7047331771567910e-03),mfarraydef(  3.0347083105033588e-04),mfarraydef(  8.3931798507141038e-06),mfarraydef(  1.6509899343992701e-07),mfarraydef(  2.2953983633288216e-09),mfarraydef(  2.2426481120489624e-11),mfarraydef(  1.5049336925781968e-13),mfarraydef(  6.2806501303798073e-16),mfarraydef(  1.2085940888554432e-18)
         };
      VReal sum1,sum2;
      const VReal x = thisVm;

      sum1 = make_float(0);
      #pragma unroll(Mhu_m)
      for (int j = Mhu_m-1; j >= 0; j--)
         sum1 = add(mfarrayuse(Mhu_a[j]),mul(x,sum1));

      sum2 = make_float(0);
      int k = Mhu_m + Mhu_l - 1;
      #pragma unroll(Mhu_l)
      for (int j = k; j >= Mhu_m; j--)
         sum2 = add(mfarrayuse(Mhu_a[j]),mul(x,sum2));
      VReal mhu = div(sum1,sum2);

      sum1 = make_float(0);
      #pragma unroll(Tau_m)
      for (int j = Tau_m-1; j >= 0; j--)
         sum1 = add(mfarrayuse(Tau_a[j]),mul(x,sum1));

      sum2 = make_float(0);
      k = Tau_m + Tau_l - 1;
      #pragma unroll(Tau_l)
      for (int j = k; j >= Tau_m; j--)
         sum2 = add(mfarrayuse(Tau_a[j]),mul(x,sum2));
      VReal tauR = div(sum1,sum2)*dt;
         ForwardEuler(hGate);
      }
   }
}
void jGate(const Real _dt, const int nSimd, const int simdBegin, const alignedRealPtr Vm, const alignedRealPtr iStim, alignedRealPtr dVm, alignedVRealPtr thisState_)
{

  const VReal dt = splat(&_dt);
   for (int ii=simdBegin; ii<nSimd+simdBegin; ++ii)
   {
      const VReal thisVm = load(&Vm[ii*SIMDOPS_FLOAT64V_WIDTH]);
      //2
      {
         const int Mhu_l = 7;
         const int Mhu_m = 9;
         const VRealArray Mhu_a[] = {
           mfarraydef(3.1263044718396443e-08),mfarraydef(  6.1534783640833612e-09),mfarraydef( -3.7646458226792461e-10),mfarraydef( -4.2451512350227289e-11),mfarraydef(  1.5216916689168596e-12),mfarraydef(  3.9843764909619204e-14),mfarraydef( -7.7509942611798927e-16),mfarraydef( -3.3646392579644030e-17),mfarraydef(  5.7566741989713017e-19),mfarraydef(  1.0000000000000000e+00),mfarraydef(  6.4169311162380532e-02),mfarraydef(  1.7346341071731517e-03),mfarraydef(  2.5211977280197058e-05),mfarraydef(  2.0770353677292379e-07),mfarraydef(  9.2161835229419180e-10),mfarraydef(  1.7368930279808849e-12)
         };
         const int Tau_l = 1;
         const int Tau_m = 13;
         const VRealArray Tau_a[] = {
           mfarraydef(5.7685269377162662e-01),mfarraydef(  3.5119791892517939e-02),mfarraydef(  9.5403324816453923e-04),mfarraydef(  1.8482121521375892e-05),mfarraydef(  2.9493061828484874e-07),mfarraydef(  2.6326436079978290e-09),mfarraydef(  9.2707225823646994e-12),mfarraydef(  3.6387910084447416e-13),mfarraydef(  8.2748066828173235e-15),mfarraydef(  2.3613021276064166e-17),mfarraydef( -8.1293559809646280e-19),mfarraydef( -7.9897103307330249e-21),mfarraydef( -2.1875746179888087e-23),mfarraydef(  1.0000000000000000e+00)
         };
      VReal sum1,sum2;
      const VReal x = thisVm;

      sum1 = make_float(0);
      #pragma unroll(Mhu_m)
      for (int j = Mhu_m-1; j >= 0; j--)
         sum1 = add(mfarrayuse(Mhu_a[j]),mul(x,sum1));

      sum2 = make_float(0);
      int k = Mhu_m + Mhu_l - 1;
      #pragma unroll(Mhu_l)
      for (int j = k; j >= Mhu_m; j--)
         sum2 = add(mfarrayuse(Mhu_a[j]),mul(x,sum2));
      VReal mhu = div(sum1,sum2);

      sum1 = make_float(0);
      #pragma unroll(Tau_m)
      for (int j = Tau_m-1; j >= 0; j--)
         sum1 = add(mfarrayuse(Tau_a[j]),mul(x,sum1));

      sum2 = make_float(0);
      k = Tau_m + Tau_l - 1;
      #pragma unroll(Tau_l)
      for (int j = k; j >= Tau_m; j--)
         sum2 = add(mfarrayuse(Tau_a[j]),mul(x,sum2));
      VReal tauR = div(sum1,sum2)*dt;
         ForwardEuler(jGate);
      }
   }
}
void Xr1Gate(const Real _dt, const int nSimd, const int simdBegin, const alignedRealPtr Vm, const alignedRealPtr iStim, alignedRealPtr dVm, alignedVRealPtr thisState_)
{

  const VReal dt = splat(&_dt);
   for (int ii=simdBegin; ii<nSimd+simdBegin; ++ii)
   {
      const VReal thisVm = load(&Vm[ii*SIMDOPS_FLOAT64V_WIDTH]);
      //3
      {
         const int Mhu_l = 5;
         const int Mhu_m = 8;
         const VRealArray Mhu_a[] = {
           mfarraydef(9.7620687511302562e-01),mfarraydef(  5.5687557065868629e-02),mfarraydef(  1.3465134489445548e-03),mfarraydef(  1.7240843920997653e-05),mfarraydef(  1.1567866757841745e-07),mfarraydef(  2.8292918822319700e-10),mfarraydef( -8.2504917317015575e-13),mfarraydef( -4.5498395933488446e-15),mfarraydef(  1.0000000000000000e+00),mfarraydef(  5.3646986850638727e-02),mfarraydef(  1.4282138385673742e-03),mfarraydef(  1.5252838752755205e-05),mfarraydef(  1.4682559892184376e-07)
         };
         const int Tau_l = 1;
         const int Tau_m = 13;
         const VRealArray Tau_a[] = {
           mfarraydef(5.4609851558432825e-03),mfarraydef(  4.3637206691412927e-04),mfarraydef(  1.9034136985586222e-05),mfarraydef(  5.4971744359155369e-07),mfarraydef(  1.2019814013692350e-08),mfarraydef(  2.0976070086137667e-10),mfarraydef(  3.0204999531382605e-12),mfarraydef(  3.6013108734580622e-14),mfarraydef(  3.9727446952762997e-16),mfarraydef(  4.3797461895814199e-18),mfarraydef(  4.2418783268539835e-20),mfarraydef(  2.6379981854029984e-22),mfarraydef(  7.6831345074599439e-25),mfarraydef(  1.0000000000000000e+00)
         };
      VReal sum1,sum2;
      const VReal x = thisVm;

      sum1 = make_float(0);
      #pragma unroll(Mhu_m)
      for (int j = Mhu_m-1; j >= 0; j--)
         sum1 = add(mfarrayuse(Mhu_a[j]),mul(x,sum1));

      sum2 = make_float(0);
      int k = Mhu_m + Mhu_l - 1;
      #pragma unroll(Mhu_l)
      for (int j = k; j >= Mhu_m; j--)
         sum2 = add(mfarrayuse(Mhu_a[j]),mul(x,sum2));
      VReal mhu = div(sum1,sum2);

      sum1 = make_float(0);
      #pragma unroll(Tau_m)
      for (int j = Tau_m-1; j >= 0; j--)
         sum1 = add(mfarrayuse(Tau_a[j]),mul(x,sum1));

      sum2 = make_float(0);
      k = Tau_m + Tau_l - 1;
      #pragma unroll(Tau_l)
      for (int j = k; j >= Tau_m; j--)
         sum2 = add(mfarrayuse(Tau_a[j]),mul(x,sum2));
      VReal tauR = div(sum1,sum2)*dt;
         ForwardEuler(Xr1Gate);

      }
   }
}
void Xr2Gate(const Real _dt, const int nSimd, const int simdBegin, const alignedRealPtr Vm, const alignedRealPtr iStim, alignedRealPtr dVm, alignedVRealPtr thisState_)
{

  const VReal dt = splat(&_dt);
   for (int ii=simdBegin; ii<nSimd+simdBegin; ++ii)
   {
      const VReal thisVm = load(&Vm[ii*SIMDOPS_FLOAT64V_WIDTH]);
      //4
      {
         const int Mhu_l = 1;
         const int Mhu_m = 10;
         const VRealArray Mhu_a[] = {
           mfarraydef(2.4919224554816769e-02),mfarraydef( -1.0132133908609619e-03),mfarraydef(  2.0115370077599340e-05),mfarraydef( -2.4787509921554170e-07),mfarraydef(  1.9157062546455413e-09),mfarraydef( -1.1856625513158431e-11),mfarraydef(  4.3373559510489725e-14),mfarraydef(  2.9410875655984132e-15),mfarraydef( -2.8535346738989432e-17),mfarraydef( -2.7319284078381429e-19),mfarraydef(  1.0000000000000000e+00)
         };
         const int Tau_l = 1;
         const int Tau_m = 10;
         const VRealArray Tau_a[] = {
           mfarraydef(3.2797890377243832e-01),mfarraydef(  1.1724387906430438e-06),mfarraydef(  3.7188870926079577e-05),mfarraydef( -3.3513027526302832e-09),mfarraydef(  7.4502277507770199e-09),mfarraydef(  1.2888652361968656e-12),mfarraydef(  8.0734422444730268e-13),mfarraydef(  8.1028416357975461e-16),mfarraydef( -1.9648667839939542e-18),mfarraydef( -3.5560477657149278e-19),mfarraydef(  1.0000000000000000e+00)
         };
      VReal sum1,sum2;
      const VReal x = thisVm;

      sum1 = make_float(0);
      #pragma unroll(Mhu_m)
      for (int j = Mhu_m-1; j >= 0; j--)
         sum1 = add(mfarrayuse(Mhu_a[j]),mul(x,sum1));

      sum2 = make_float(0);
      int k = Mhu_m + Mhu_l - 1;
      #pragma unroll(Mhu_l)
      for (int j = k; j >= Mhu_m; j--)
         sum2 = add(mfarrayuse(Mhu_a[j]),mul(x,sum2));
      VReal mhu = div(sum1,sum2);

      sum1 = make_float(0);
      #pragma unroll(Tau_m)
      for (int j = Tau_m-1; j >= 0; j--)
         sum1 = add(mfarrayuse(Tau_a[j]),mul(x,sum1));

      sum2 = make_float(0);
      k = Tau_m + Tau_l - 1;
      #pragma unroll(Tau_l)
      for (int j = k; j >= Tau_m; j--)
         sum2 = add(mfarrayuse(Tau_a[j]),mul(x,sum2));
      VReal tauR = div(sum1,sum2)*dt;
         ForwardEuler(Xr2Gate);

      }
   }
}
void XsGate(const Real _dt, const int nSimd, const int simdBegin, const alignedRealPtr Vm, const alignedRealPtr iStim, alignedRealPtr dVm, alignedVRealPtr thisState_)
{

  const VReal dt = splat(&_dt);
   for (int ii=simdBegin; ii<nSimd+simdBegin; ++ii)
   {
      const VReal thisVm = load(&Vm[ii*SIMDOPS_FLOAT64V_WIDTH]);
      //5
      {
         const int Mhu_l = 5;
         const int Mhu_m = 5;
         const VRealArray Mhu_a[] = {
           mfarraydef(5.8834144047581471e-01),mfarraydef(  1.8554918971499884e-02),mfarraydef(  2.3769494810239468e-04),mfarraydef(  1.4558400579866515e-06),mfarraydef(  3.5608483256093079e-09),mfarraydef(  1.0000000000000000e+00),mfarraydef(  2.1263634327928394e-03),mfarraydef(  5.2685461356465763e-04),mfarraydef( -1.2603117542017038e-06),mfarraydef(  1.5082411859988069e-08)
         };
         const int Tau_l = 6;
         const int Tau_m = 9;
         const VRealArray Tau_a[] = {
           mfarraydef(1.2782667970990670e-03),mfarraydef( -5.3480687790940232e-05),mfarraydef(  4.1463175377037539e-06),mfarraydef( -7.9427267607107334e-08),mfarraydef(  2.3181534592617566e-09),mfarraydef( -1.8597152793072996e-11),mfarraydef(  3.4109130605190490e-13),mfarraydef(  2.0323839081334259e-15),mfarraydef(  4.3130861032059193e-18),mfarraydef(  1.0000000000000000e+00),mfarraydef(  4.9865573111495263e-03),mfarraydef(  1.0343375562869115e-03),mfarraydef(  1.3748473468950345e-05),mfarraydef(  3.3217061243887820e-07),mfarraydef( -2.7950790321181472e-09)
         };
      VReal sum1,sum2;
      const VReal x = thisVm;

      sum1 = make_float(0);
      #pragma unroll(Mhu_m)
      for (int j = Mhu_m-1; j >= 0; j--)
         sum1 = add(mfarrayuse(Mhu_a[j]),mul(x,sum1));

      sum2 = make_float(0);
      int k = Mhu_m + Mhu_l - 1;
      #pragma unroll(Mhu_l)
      for (int j = k; j >= Mhu_m; j--)
         sum2 = add(mfarrayuse(Mhu_a[j]),mul(x,sum2));
      VReal mhu = div(sum1,sum2);

      sum1 = make_float(0);
      #pragma unroll(Tau_m)
      for (int j = Tau_m-1; j >= 0; j--)
         sum1 = add(mfarrayuse(Tau_a[j]),mul(x,sum1));

      sum2 = make_float(0);
      k = Tau_m + Tau_l - 1;
      #pragma unroll(Tau_l)
      for (int j = k; j >= Tau_m; j--)
         sum2 = add(mfarrayuse(Tau_a[j]),mul(x,sum2));
      VReal tauR = div(sum1,sum2)*dt;
         ForwardEuler(XsGate);
      }
   }
}
void rGate(const Real _dt, const int nSimd, const int simdBegin, const alignedRealPtr Vm, const alignedRealPtr iStim, alignedRealPtr dVm, alignedVRealPtr thisState_)
{

  const VReal dt = splat(&_dt);
   for (int ii=simdBegin; ii<nSimd+simdBegin; ++ii)
   {
      const VReal thisVm = load(&Vm[ii*SIMDOPS_FLOAT64V_WIDTH]);
      //6
      {
         const int Mhu_l = 5;
         const int Mhu_m = 8;
         const VRealArray Mhu_a[] = {
           mfarraydef(3.4450952356981063e-02),mfarraydef(  3.2569778153723060e-03),mfarraydef(  1.4015137277741737e-04),mfarraydef(  3.5177311829244476e-06),mfarraydef(  5.4744904248961321e-08),mfarraydef(  5.1978297648799841e-10),mfarraydef(  2.7509401026136440e-12),mfarraydef(  6.2017518287576868e-15),mfarraydef(  1.0000000000000000e+00),mfarraydef( -6.6423744780337163e-02),mfarraydef(  2.2756546586311602e-03),mfarraydef( -3.0510539067819270e-05),mfarraydef(  3.1936830897476799e-07)
         };
         const int Tau_l = 5;
         const int Tau_m = 7;
         const VRealArray Tau_a[] = {
           mfarraydef(2.1250755813646321e-01),mfarraydef(  5.4375121806696314e-03),mfarraydef(  1.2855334968507940e-04),mfarraydef(  1.8435064622728614e-06),mfarraydef(  2.2201936942212559e-08),mfarraydef(  1.6208374761956396e-10),mfarraydef(  7.0475957646256506e-13),mfarraydef(  1.0000000000000000e+00),mfarraydef( -1.1298206006665951e-02),mfarraydef(  1.9108836371246739e-05),mfarraydef(  3.8822652158898901e-06),mfarraydef(  2.4266538879422829e-08)
         };
      VReal sum1,sum2;
      const VReal x = thisVm;

      sum1 = make_float(0);
      #pragma unroll(Mhu_m)
      for (int j = Mhu_m-1; j >= 0; j--)
         sum1 = add(mfarrayuse(Mhu_a[j]),mul(x,sum1));

      sum2 = make_float(0);
      int k = Mhu_m + Mhu_l - 1;
      #pragma unroll(Mhu_l)
      for (int j = k; j >= Mhu_m; j--)
         sum2 = add(mfarrayuse(Mhu_a[j]),mul(x,sum2));
      VReal mhu = div(sum1,sum2);

      sum1 = make_float(0);
      #pragma unroll(Tau_m)
      for (int j = Tau_m-1; j >= 0; j--)
         sum1 = add(mfarrayuse(Tau_a[j]),mul(x,sum1));

      sum2 = make_float(0);
      k = Tau_m + Tau_l - 1;
      #pragma unroll(Tau_l)
      for (int j = k; j >= Tau_m; j--)
         sum2 = add(mfarrayuse(Tau_a[j]),mul(x,sum2));
      VReal tauR = div(sum1,sum2)*dt;
         ForwardEuler(rGate);
      }
   }
}
void dGate(const Real _dt, const int nSimd, const int simdBegin, const alignedRealPtr Vm, const alignedRealPtr iStim, alignedRealPtr dVm, alignedVRealPtr thisState_)
{

  const VReal dt = splat(&_dt);
   for (int ii=simdBegin; ii<nSimd+simdBegin; ++ii)
   {
      const VReal thisVm = load(&Vm[ii*SIMDOPS_FLOAT64V_WIDTH]);
      //7
      {
         const int Mhu_l = 5;
         const int Mhu_m = 7;
         const VRealArray Mhu_a[] = {
           mfarraydef(            7.4398777641576486e-01),mfarraydef(  4.4326834369266523e-02),mfarraydef(  1.1359075050132548e-03),mfarraydef(  1.5872707499276087e-05),mfarraydef(  1.2622142062210935e-07),mfarraydef(  5.3559505029241141e-10),mfarraydef(  9.3621878808126576e-13),mfarraydef(  1.0000000000000000e+00),mfarraydef(  2.5428070655441278e-02),mfarraydef(  1.7701065218545982e-03),mfarraydef(  3.7348639823389864e-06),mfarraydef(  2.5522660570009651e-07)
         };
         const int Tau_l = 7;
         const int Tau_m = 10;
         const VRealArray Tau_a[] = {
           mfarraydef(1.5068318180921829e+00),mfarraydef(  7.8993621885004611e-02),mfarraydef(  6.9687982613764387e-03),mfarraydef(  1.1019128271535092e-04),mfarraydef(  2.3993512027253813e-06),mfarraydef(  1.1356515133023819e-08),mfarraydef(  3.1836750643625883e-10),mfarraydef( -1.6718952430054101e-12),mfarraydef(  2.5579267994953318e-14),mfarraydef(  1.0819050682737564e-16),mfarraydef(  1.0000000000000000e+00),mfarraydef( -6.8101965713816826e-02),mfarraydef(  3.7430651166330204e-03),mfarraydef(  6.6379282851072469e-07),mfarraydef(  9.3098141480209233e-07),mfarraydef(  3.2614589730112466e-08),mfarraydef(  4.4763500756297451e-10)
         };
      VReal sum1,sum2;
      const VReal x = thisVm;

      sum1 = make_float(0);
      #pragma unroll(Mhu_m)
      for (int j = Mhu_m-1; j >= 0; j--)
         sum1 = add(mfarrayuse(Mhu_a[j]),mul(x,sum1));

      sum2 = make_float(0);
      int k = Mhu_m + Mhu_l - 1;
      #pragma unroll(Mhu_l)
      for (int j = k; j >= Mhu_m; j--)
         sum2 = add(mfarrayuse(Mhu_a[j]),mul(x,sum2));
      VReal mhu = div(sum1,sum2);

      sum1 = make_float(0);
      #pragma unroll(Tau_m)
      for (int j = Tau_m-1; j >= 0; j--)
         sum1 = add(mfarrayuse(Tau_a[j]),mul(x,sum1));

      sum2 = make_float(0);
      k = Tau_m + Tau_l - 1;
      #pragma unroll(Tau_l)
      for (int j = k; j >= Tau_m; j--)
         sum2 = add(mfarrayuse(Tau_a[j]),mul(x,sum2));
      VReal tauR = div(sum1,sum2)*dt;
         ForwardEuler(dGate);
      }
   }
}
void fGate(const Real _dt, const int nSimd, const int simdBegin, const alignedRealPtr Vm, const alignedRealPtr iStim, alignedRealPtr dVm, alignedVRealPtr thisState_)
{

  const VReal dt = splat(&_dt);
   for (int ii=simdBegin; ii<nSimd+simdBegin; ++ii)
   {
      const VReal thisVm = load(&Vm[ii*SIMDOPS_FLOAT64V_WIDTH]);
      //8
      {
         const int Mhu_l = 5;
         const int Mhu_m = 8;
         const VRealArray Mhu_a[] = {
           mfarraydef(5.4304780098652715e-02),mfarraydef( -4.4962960006268609e-03),mfarraydef(  1.7191393153273102e-04),mfarraydef( -3.9003819270270578e-06),mfarraydef(  5.4530196576796994e-08),mfarraydef( -4.1503002319748868e-10),mfarraydef(  8.3481595187685469e-13),mfarraydef(  5.8904944849510583e-15),mfarraydef(  1.0000000000000000e+00),mfarraydef(  5.2318145763546209e-02),mfarraydef(  1.6278824780069224e-03),mfarraydef(  1.6010969045398731e-05),mfarraydef(  1.9936720027259198e-07)
         };
         const int Tau_l = 13;
         const int Tau_m = 11;
         const VRealArray Tau_a[] = {
           mfarraydef(8.7301971950165709e-03),mfarraydef(  1.0619909606051185e-03),mfarraydef(  7.4439522770957243e-05),mfarraydef(  3.5189522703546544e-06),mfarraydef(  1.2763275154107722e-07),mfarraydef(  3.9329666177349803e-09),mfarraydef(  1.0963505894631244e-10),mfarraydef(  2.5405559223755679e-12),mfarraydef(  4.1926387533661996e-14),mfarraydef(  4.1426318225610267e-16),mfarraydef(  1.8370569627175333e-18),mfarraydef(  1.0000000000000000e+00),mfarraydef(  5.3496360858059103e-02),mfarraydef(  1.0576677271235941e-02),mfarraydef(  6.0292981147301659e-04),mfarraydef(  2.4470545201049391e-05),mfarraydef(  8.8342556285406000e-07),mfarraydef(  2.5756751722220191e-08),mfarraydef(  5.6840520101297543e-10),mfarraydef(  8.9563864907355775e-12),mfarraydef(  8.8109045583454404e-14),mfarraydef(  4.1790988048136637e-16),mfarraydef(  3.0481930918230229e-19),mfarraydef(  7.4195798494855779e-22)
         };
      VReal sum1,sum2;
      const VReal x = thisVm;

      sum1 = make_float(0);
      #pragma unroll(Mhu_m)
      for (int j = Mhu_m-1; j >= 0; j--)
         sum1 = add(mfarrayuse(Mhu_a[j]),mul(x,sum1));

      sum2 = make_float(0);
      int k = Mhu_m + Mhu_l - 1;
      #pragma unroll(Mhu_l)
      for (int j = k; j >= Mhu_m; j--)
         sum2 = add(mfarrayuse(Mhu_a[j]),mul(x,sum2));
      VReal mhu = div(sum1,sum2);

      sum1 = make_float(0);
      #pragma unroll(Tau_m)
      for (int j = Tau_m-1; j >= 0; j--)
         sum1 = add(mfarrayuse(Tau_a[j]),mul(x,sum1));

      sum2 = make_float(0);
      k = Tau_m + Tau_l - 1;
      #pragma unroll(Tau_l)
      for (int j = k; j >= Tau_m; j--)
         sum2 = add(mfarrayuse(Tau_a[j]),mul(x,sum2));
      VReal tauR = div(sum1,sum2)*dt;
         ForwardEuler(fGate);
      }
   }
}
void f2Gate(const Real _dt, const int nSimd, const int simdBegin, const alignedRealPtr Vm, const alignedRealPtr iStim, alignedRealPtr dVm, alignedVRealPtr thisState_)
{

  const VReal dt = splat(&_dt);
   for (int ii=simdBegin; ii<nSimd+simdBegin; ++ii)
   {
      const VReal thisVm = load(&Vm[ii*SIMDOPS_FLOAT64V_WIDTH]);
      //9
      {
         const int Mhu_l = 8;
         const int Mhu_m = 5;
         const VRealArray Mhu_a[] = {
           mfarraydef(3.3448475117101473e-01),mfarraydef(  1.7908684328425843e-02),mfarraydef(  4.6473146231373311e-04),mfarraydef(  4.9947394623057196e-06),mfarraydef(  4.6013519663601120e-08),mfarraydef(  1.0000000000000000e+00),mfarraydef(  5.5442666103601546e-02),mfarraydef(  1.3608399868819414e-03),mfarraydef(  1.6306821023781025e-05),mfarraydef(  1.2063672477786111e-07),mfarraydef(  1.7672562299303354e-10),mfarraydef( -5.4528541905004782e-13),mfarraydef( -2.8928059365721159e-15)
         };
         const int Tau_l = 11;
         const int Tau_m = 12;
         const VRealArray Tau_a[] = {
           mfarraydef(3.0215260344315119e-02),mfarraydef(  2.9905565403315244e-03),mfarraydef(  1.8698675038604531e-04),mfarraydef(  7.8302945745729333e-06),mfarraydef(  2.4353851952274140e-07),mfarraydef(  5.9716586064763031e-09),mfarraydef(  1.4320157611618665e-10),mfarraydef(  3.5798223309278773e-12),mfarraydef(  7.3139244140308910e-14),mfarraydef(  9.1797821675112518e-16),mfarraydef(  5.5127688299509001e-18),mfarraydef(  6.0160642283988931e-21),mfarraydef(  1.0000000000000000e+00),mfarraydef( -8.8631007579573415e-02),mfarraydef(  5.6222208830099462e-03),mfarraydef(  9.0331117197141856e-05),mfarraydef(  8.5481830250918979e-07),mfarraydef(  2.2780719183249479e-07),mfarraydef(  7.0378379382184305e-09),mfarraydef(  1.0211239994063565e-10),mfarraydef(  1.8042686984162349e-12),mfarraydef(  2.9931381582699879e-14),mfarraydef(  2.1128053668752277e-16)
         };
      VReal sum1,sum2;
      const VReal x = thisVm;

      sum1 = make_float(0);
      #pragma unroll(Mhu_m)
      for (int j = Mhu_m-1; j >= 0; j--)
         sum1 = add(mfarrayuse(Mhu_a[j]),mul(x,sum1));

      sum2 = make_float(0);
      int k = Mhu_m + Mhu_l - 1;
      #pragma unroll(Mhu_l)
      for (int j = k; j >= Mhu_m; j--)
         sum2 = add(mfarrayuse(Mhu_a[j]),mul(x,sum2));
      VReal mhu = div(sum1,sum2);

      sum1 = make_float(0);
      #pragma unroll(Tau_m)
      for (int j = Tau_m-1; j >= 0; j--)
         sum1 = add(mfarrayuse(Tau_a[j]),mul(x,sum1));

      sum2 = make_float(0);
      k = Tau_m + Tau_l - 1;
      #pragma unroll(Tau_l)
      for (int j = k; j >= Tau_m; j--)
         sum2 = add(mfarrayuse(Tau_a[j]),mul(x,sum2));
      VReal tauR = div(sum1,sum2)*dt;
         ForwardEuler(f2Gate);
      }
   }
}
void jLGate(const Real _dt, const int nSimd, const int simdBegin, const alignedRealPtr Vm, const alignedRealPtr iStim, alignedRealPtr dVm, alignedVRealPtr thisState_)
{

  const VReal dt = splat(&_dt);
   for (int ii=simdBegin; ii<nSimd+simdBegin; ++ii)
   {
      const VReal thisVm = load(&Vm[ii*SIMDOPS_FLOAT64V_WIDTH]);
      //10
      {
         const int Mhu_l = 6;
         const int Mhu_m = 12;
         const VRealArray Mhu_a[] = {
           mfarraydef(3.6932725836170428e-09),mfarraydef(  2.5736169972850485e-10),mfarraydef( -7.5177270851662063e-11),mfarraydef( -1.7783383158540549e-12),mfarraydef(  2.3780959584764758e-13),mfarraydef(  3.9995914648833690e-15),mfarraydef( -2.6002507057823451e-16),mfarraydef( -4.0173750386716292e-18),mfarraydef(  1.1567747025479995e-19),mfarraydef(  1.8190095403795220e-21),mfarraydef( -1.7772014587013272e-23),mfarraydef( -3.0513026345262148e-25),mfarraydef(  1.0000000000000000e+00),mfarraydef(  4.7761481528554264e-02),mfarraydef(  9.1744037118706855e-04),mfarraydef(  8.8560239615590481e-06),mfarraydef(  4.2943755440989806e-08),mfarraydef(  8.3656852902593096e-11)
         };
         const int Tau_l = 1;
         const int Tau_m = 1;
         const VRealArray Tau_a[] = {
           mfarraydef(1.4925373134328358e-03),mfarraydef(  1.0000000000000000e+00)
         };
      VReal sum1,sum2;
      const VReal x = thisVm;

      sum1 = make_float(0);
      #pragma unroll(Mhu_m)
      for (int j = Mhu_m-1; j >= 0; j--)
         sum1 = add(mfarrayuse(Mhu_a[j]),mul(x,sum1));

      sum2 = make_float(0);
      int k = Mhu_m + Mhu_l - 1;
      #pragma unroll(Mhu_l)
      for (int j = k; j >= Mhu_m; j--)
         sum2 = add(mfarrayuse(Mhu_a[j]),mul(x,sum2));
      VReal mhu = div(sum1,sum2);

      sum1 = make_float(0);
      #pragma unroll(Tau_m)
      for (int j = Tau_m-1; j >= 0; j--)
         sum1 = add(mfarrayuse(Tau_a[j]),mul(x,sum1));

      sum2 = make_float(0);
      k = Tau_m + Tau_l - 1;
      #pragma unroll(Tau_l)
      for (int j = k; j >= Tau_m; j--)
         sum2 = add(mfarrayuse(Tau_a[j]),mul(x,sum2));
      VReal tauR = div(sum1,sum2)*dt;
         ForwardEuler(jLGate);
      }
   }
}
void sGate(const Real _dt, const int nSimd, const int simdBegin, const alignedRealPtr Vm, const alignedRealPtr iStim, alignedRealPtr dVm, alignedVRealPtr thisState_)
{

  const VReal dt = splat(&_dt);
   for (int ii=simdBegin; ii<nSimd+simdBegin; ++ii)
   {
      const VReal thisVm = load(&Vm[ii*SIMDOPS_FLOAT64V_WIDTH]);
      //11
      {
         const int Mhu_l = 7;
         const int Mhu_m = 8;
         const VRealArray Mhu_a[] = {
           mfarraydef(3.6866939888900814e-03),mfarraydef( -4.5574539113452475e-04),mfarraydef(  2.6788146214637530e-05),mfarraydef( -9.9102133439606097e-07),mfarraydef(  2.5613910193807883e-08),mfarraydef( -4.6400777680998212e-10),mfarraydef(  5.2352368345974410e-12),mfarraydef( -2.6731132916017718e-14),mfarraydef(  1.0000000000000000e+00),mfarraydef(  7.5601232522433295e-02),mfarraydef(  2.5744765696710309e-03),mfarraydef(  4.7389711229137498e-05),mfarraydef(  5.8624244628409614e-07),mfarraydef(  3.4873125981144344e-09),mfarraydef(  2.0871941515527373e-11)
         };
         const int Tau_l = 6;
         const int Tau_m = 7;
         const VRealArray Tau_a[] = {
           mfarraydef(5.2008879668973697e-02),mfarraydef(  2.9402318950257256e-03),mfarraydef(  9.6193142778593343e-05),mfarraydef(  1.9694354926101320e-06),mfarraydef(  2.5766968301605511e-08),mfarraydef(  1.9309877695248136e-10),mfarraydef(  6.5471492087484267e-13),mfarraydef(  1.0000000000000000e+00),mfarraydef( -2.1715782138805460e-02),mfarraydef(  2.0915626226951107e-03),mfarraydef(  1.8073629807905746e-07),mfarraydef(  2.1067307291025196e-07),mfarraydef(  2.7434857789930993e-09)
         };
      VReal sum1,sum2;
      const VReal x = thisVm;

      sum1 = make_float(0);
      #pragma unroll(Mhu_m)
      for (int j = Mhu_m-1; j >= 0; j--)
         sum1 = add(mfarrayuse(Mhu_a[j]),mul(x,sum1));

      sum2 = make_float(0);
      int k = Mhu_m + Mhu_l - 1;
      #pragma unroll(Mhu_l)
      for (int j = k; j >= Mhu_m; j--)
         sum2 = add(mfarrayuse(Mhu_a[j]),mul(x,sum2));
      VReal mhu = div(sum1,sum2);

      sum1 = make_float(0);
      #pragma unroll(Tau_m)
      for (int j = Tau_m-1; j >= 0; j--)
         sum1 = add(mfarrayuse(Tau_a[j]),mul(x,sum1));

      sum2 = make_float(0);
      k = Tau_m + Tau_l - 1;
      #pragma unroll(Tau_l)
      for (int j = k; j >= Tau_m; j--)
         sum2 = add(mfarrayuse(Tau_a[j]),mul(x,sum2));
      VReal tauR = div(sum1,sum2)*dt;
         ForwardEuler(sGate);
      }
   }
}

 typedef void (*gateFuncType)(const Real _dt, const int nSimd, const int simdBegin, const alignedRealPtr Vm, const alignedRealPtr iStim, alignedRealPtr dVm, alignedVRealPtr thisState_);


alignedVRealPtr initState(const double initCond, const int nsize) {
  alignedVRealPtr retval = new VReal[nsize];
  for (int ii=0; ii<nsize; ii++) {
    native_vector_type tmp = make_float(initCond);
    retval[ii] = tmp;
//    retval[ii] = make_float(initCond);
  }
  return retval;
}

int main(void) {

   int nCells = (1 << 8); //approx 1M cells.
   int numThreads = omp_get_max_threads();
   int ncpu=omp_get_num_procs();

   nCells -= (nCells % (numThreads*SIMDOPS_FLOAT64V_WIDTH));

   cout << numThreads << endl;
   cout << ncpu << endl;
   cout << nCells << endl;

   int numGates = 12;
   int threadsPerGate[numGates];
   for (int ii=0; ii<numGates; ii++) {
     threadsPerGate[ii] = numThreads/numGates;
     threadsPerGate[ii] += (ii < (numThreads % numGates));
   }
   int threadExtents[numGates+1];
   threadExtents[0] = 0;
   for (int ii=1; ii<=numGates; ii++) {
     threadExtents[ii] = threadsPerGate[ii-1]+threadExtents[ii-1];
   }
   gateFuncType gateFunc[numGates];
   gateFunc[ 0] = f2Gate;
   gateFunc[ 1] = fGate;
   gateFunc[ 2] = dGate;
   gateFunc[ 3] = mGate;
   gateFunc[ 4] = jGate;
   gateFunc[ 5] = hGate;
   gateFunc[ 6] = rGate;
   gateFunc[ 7] = sGate;
   gateFunc[ 8] = Xr1Gate;
   gateFunc[ 9] = Xr2Gate;
   gateFunc[10] = XsGate;
   gateFunc[11] = jLGate;


   vectorSimd(Real) Vm(nCells);
   vectorSimd(Real) dVm(nCells);
   vectorSimd(Real) iStim(nCells);
   int nState = (nCells+SIMDOPS_FLOAT64V_WIDTH-1)/SIMDOPS_FLOAT64V_WIDTH;
   State state;
   cout << SIMDOPS_FLOAT64V_WIDTH <<endl;

   {
      const Real initVm = -86.709;
      Vm.assign(nCells, initVm);
      dVm.assign(nCells, 0);

      const Real pcnst_2 = 96485.3415;
      const Real pcnst_3 = 0.185;
      const Real pcnst_4 = 0.016404;
      const Real c9 = -pcnst_3/(pcnst_4*pcnst_2);
      const Real K_i = 138.4;
      state.dVK_i    =initState(K_i/c9+initVm, nState);
      state.Na_i     =initState(10.355  , nState);
      state.Ca_i     =initState(0.00013 , nState);
      state.Ca_ss    =initState(0.00036 , nState);
      state.Ca_sr    =initState(3.715   , nState);
      state.R_prime  =initState(0.9068  , nState);
      state.fCass    =initState(0.9953  , nState);
      state.Xr1Gate  =initState(0.00448 , nState);
      state.Xr2Gate  =initState(0.476   , nState);
      state.XsGate   =initState(0.0087  , nState);
      state.mGate    =initState(0.00155 , nState);
      state.hGate    =initState(0.7573  , nState);
      state.jGate    =initState(0.7225  , nState);
      state.rGate    =initState(2.235e-8, nState);
      state.dGate    =initState(3.164e-5, nState);
      state.fGate    =initState(0.8009  , nState);
      state.f2Gate   =initState(0.9778  , nState);
      state.sGate    =initState(0.3212  , nState);
      state.jLGate   =initState(0.066   , nState);
   }
   alignedVRealPtr gateArray[numGates];
   gateArray[ 0] = state.f2Gate;
   gateArray[ 1] = state.fGate;
   gateArray[ 2] = state.dGate;
   gateArray[ 3] = state.mGate;
   gateArray[ 4] = state.jGate;
   gateArray[ 5] = state.hGate;
   gateArray[ 6] = state.rGate;
   gateArray[ 7] = state.sGate;
   gateArray[ 8] = state.Xr1Gate;
   gateArray[ 9] = state.Xr2Gate;
   gateArray[10] = state.XsGate;
   gateArray[11] = state.jLGate;

   HPM_Init();
   HPM_Start("main");

   Real dt=0.02;
   int nSteps = 10000;
#pragma omp parallel
   {
   int tid = omp_get_thread_num();
   int procid = tid;
   int coreid = tid%ncpu;
   int hwthreadid = tid/ncpu;
   
   int myNCells = nCells/numThreads;

   int myCellsBegin = (nCells/numThreads) * tid;
   int myCellsEnd = myCellsBegin+myNCells;

   int myNState = nState/numThreads;
   int myStateBegin = myNState * tid;

   int myGate;
   for (int ii=0; ii<numGates; ii++) {
     if (tid < threadExtents[ii+1]) {
       myGate=ii;
       break;
     }
   }
   int tidMyGate = tid - threadExtents[myGate];
   int threadsPerMyGate = threadsPerGate[myGate];

   int nSimd = nCells/SIMDOPS_FLOAT64V_WIDTH;
   int myGateNSimd = nSimd/threadsPerMyGate + (tidMyGate < (nSimd % threadsPerMyGate)?
                                               1 : 0);
   int myGateSimdBegin = (nSimd/threadsPerMyGate) * tidMyGate
     +
     min(tidMyGate, nSimd % threadsPerMyGate)
     ;
   int myGateNCells = myGateNSimd * SIMDOPS_FLOAT64V_WIDTH;
   int myGateCellsBegin = myGateSimdBegin * SIMDOPS_FLOAT64V_WIDTH;
   
   int myGateStateBegin = myGateSimdBegin;


   /*
   #pragma omp critical
   {
      cout << tid 
           << "  " << myNCells 
           << "  " << myCellsBegin
           << "  " << myCellsEnd
           << "  " << &state[0]+myCellsBegin
           << "  " << &Vm[0]+myCellsBegin
           << "  " << &dVm[0]+myCellsBegin
           << "  " << &iStim[0]+myCellsBegin
           <<endl;
   }
   #pragma omp barrier
   */

   cout.precision(16);
   for (int ii=0; ii<nSteps+1; ii++) {
      if (tid==0 && (ii % 2000) == 0) {
         cout << dt*ii << "\t"
              << Vm[1] << "\t"
              << toDouble(state.f2Gate[1]) << "\t"
              << toDouble(state.fGate[1]) << "\t"
              << toDouble(state.dGate[1]) << "\t"
              << toDouble(state.mGate[1]) << "\t"
              << toDouble(state.jGate[1]) << "\t"
              << toDouble(state.hGate[1]) << "\t"
              << toDouble(state.rGate[1]) << "\t"
              << toDouble(state.sGate[1]) << "\t"
              << toDouble(state.Xr1Gate[1]) << "\t"
              << toDouble(state.Xr2Gate[1]) << "\t"
              << toDouble(state.XsGate[1]) << "\t"
              << toDouble(state.jLGate[1]) << "\t"
              << toDouble(state.Na_i[1]) << "\t"
              << toDouble(state.Ca_i[1]) << "\t"
              << toDouble(state.Ca_ss[1]) << "\t"
              << toDouble(state.Ca_sr[1]) << "\t"
              << toDouble(state.fCass[1]) << "\t"
              << toDouble(state.dVK_i[1]) << "\t"
              << toDouble(state.R_prime[1]) << endl;
      }
      if (ii==nSteps) { break; }

      Real istimVal = 0;
      if (ii < 50) {
         istimVal = -60;
      }
      for (int icell=myCellsBegin; icell < myCellsEnd; icell++) {
         iStim[icell] = istimVal;
      }
      nonGates(dt, myNState, myStateBegin,
               (&Vm[0]), 
               (&iStim[0]), 
               (&dVm[0]), 
               &state
               );
      #pragma omp barrier


      gateFunc[myGate](dt, myGateNSimd, myGateSimdBegin, 
                       (&Vm[0]), 
                       (&iStim[0]), 
                       (&dVm[0]), 
                       gateArray[myGate]);

      #pragma omp barrier

      for (int icell=myCellsBegin; icell<myCellsEnd; icell++) {
         Vm[icell] += dt*(dVm[icell]-istimVal);
      }
   }

   }


   HPM_Stop("main");
   cout << "print hpm output" << endl;
   HPM_Print();

   cout << globalNonGate << endl;
   cout << globalGate << endl;

   return 0;
}
