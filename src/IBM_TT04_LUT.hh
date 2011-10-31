#ifndef IBM_TT04_LUT_HH
#define IBM_TT04_LUT_HH

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;


#define STIMULATION
#ifdef STIMULATION
const int RangeTab=5000;
const int DivisionTab=10;
#else
const int RangeTab=400;
const int DivisionTab=10;
#endif

const int V_RESOLUTION=RangeTab*DivisionTab;
const int RangeTabhalf=RangeTab/2;
const double dDivisionTab=1.0/DivisionTab;


#define TT04_STATE_VARIABLES 17

enum TT04StateVariables
{
  tt04_V,
  tt04_Nai,
  tt04_Ki,
  tt04_Cai,
  tt04_xr1,
  tt04_xr2,
  tt04_xs,
  tt04_m,
  tt04_h,
  tt04_j,
  tt04_d,
  tt04_f,
  tt04_fCa,
  tt04_s,
  tt04_r,
  tt04_CaSR,
  tt04_g
};


typedef struct
{

  double R; //            // gas constant                     - J*K^(-1)*mol^(-1)
  double T; //             // temperature                      - K
  double F; //           // Faraday constant                 - C/mmol
  double Cm; //                // Cell capacitance per unit surface area  - microF/cm^2
  double VC; //            // Cytoplasmic volume               - microm^3
  double pKNa; // 0.03          // relative IKs permeability to Na+
  double Nao; // 140            // Extracellular Na+ concentration  - mM
  double Ko; // 5.4             // Extracellular K+ concentration   - mM
  double Cao; // 2              // Extracellular Ca2+ concentration - mM
  double GK1; // 5.405          // Maximal IK1 conductance          - nS/pF
  double GKr; // 0.096          // Maximal IKr conductance          - nS/pF
  double GKs_epi_endo; // 0.245 // Maximal IKs (epi/endo) conductance - nS/pF
  double GKs_M; // 0.062        // Maximal IKs (M) conductance        - nS/pF
  double GKs;
  double GNa; // 14.838         // Maximal INa conductance          - nS/pF
  double GbNa; // 0.00029       // Maximal IbNa conductance               - nS/pF
  double GCaL; // 0.000175      // = 1.75*10^(-4) != {0.10662224 = 1.75^(-4)} Maximal ICaL conductance - cm^3*microF^(-1)*s^(-1)
  double GbCa; // 0.000592      // Maximal IbCa conductance               - nS/pF
  double Gto_epi_M; // 0.294    // Maximal Gto (epi/M) conductance  - nS/pF
  double Gto_endo; // 0.073     // Maximal Gto (endo) conductance   - nS/pF
  double Gto;
  double PNaK; // 1.362         // Maximal INaK                           - pA/pF
  double KmK; // 1              // Ko half-saturation constant of INaK    - mM
  double KmNa; // 40            // Nai half-saturation constant of INaK   - mM
  double kNaCa; // 1000         // Maximal INaCa                          - pA/pF
  double ksat; // 0.1           // Saturation factor for INaCa
  double alpha; // 2.5          // Factor enhancing outward nature of INa CA
  double gamma; // 0.35         // Voltage dependence parameter of INaCa 
  double KmCa; // 1.38          // Cai half saturation constant for INaCa - mM
  double KmNai; // 87.5         // Nai half saturation constant for INaCa - mM
  double GpCa; // 0.025         // Maximal IpCa conductance               - nS/pF    IBT g_pCa .825
  double KpCa; // 0.0005        // Cai half-saturation constant of IpCa   - mM
  double GpK; // 0.0146         // Maximal IpK conductance                - nS/pF
  double arel; // 16.464        // Maximal CaSR-dependent Irel            - mM/s
  double brel; // 0.25          // CaSR half-saturation constant of Irel  - mM
  double crel; // 8.232         // Maximal CaSR-independent Irel          - mM/s
  double Kup; // 0.00025        // Half-saturation constant of Iup        - mM
  double Vleak; // 0.00008      // Maximal Ileak                          - ms^(-1)
  double Vmaxup; // 0.000425    // Maximal Iup                            - mM/ms
  double Bufc; // 0.15          // Total cytoplasmic buffer concentration - mM
  double Kbufc; // 0.001        // Cai half-saturation constant for cytoplasmic buffer concentraton - mM
  double Bufsr; // 10           // Total sarcoplasmic buffer concentration - mM
  double Kbufsr; // 0.3         // CaSR half-saturation constant for sarcoplasmic buffer - mM
  double VSR; // 1094           // Sarcoplasmic reticulum volume    - microm^3


  double RToverF;
  double inverseRToverF;
  double KopKNaNao;
  double CaL_4FF_RT;
  double inverseVcF;
  double inverse2VcF;
  double KmNai3;
  double Nao3;


/*
  //not used in single cell computation
  double S; // 0.2              // surface to volume ratio          - microm^(-1)
  double rho; // 162            // cellular resistivity             - Ohm*cm
*/

} TT04ParameterStruct;


class IBM_TT04_LUT
{
  public:
  
     // IBM parameter struct
     static char TT04LUT_EVfilename[255];// = "/bgpusr1/mr198/IBT/src/CellModel/ev_forceParameter/Flaim0910MutA/288_100.ev";
     static TT04ParameterStruct TT04LUT_ttps; //fps stands for TTParameterStruct

     static double TT04LUT_Xr1_inf[V_RESOLUTION];
     static double TT04LUT_tau_Xr1[V_RESOLUTION];
     static double TT04LUT_Xr2_inf[V_RESOLUTION];
     static double TT04LUT_tau_Xr2[V_RESOLUTION];
     static double TT04LUT_Xs_inf[V_RESOLUTION];
     static double TT04LUT_tau_Xs[V_RESOLUTION];
     static double TT04LUT_m_inf[V_RESOLUTION];
     static double TT04LUT_tau_m[V_RESOLUTION];
     static double TT04LUT_h_inf[V_RESOLUTION];
     static double TT04LUT_tau_h[V_RESOLUTION];
     static double TT04LUT_j_inf[V_RESOLUTION];
     static double TT04LUT_tau_j[V_RESOLUTION];
     static double TT04LUT_d_inf[V_RESOLUTION];
     static double TT04LUT_tau_d[V_RESOLUTION];
     static double TT04LUT_f_inf[V_RESOLUTION];
     static double TT04LUT_tau_f[V_RESOLUTION];
     static double TT04LUT_CaL_P1[V_RESOLUTION];
     static double TT04LUT_CaL_P2[V_RESOLUTION];
     static double TT04LUT_EpiM_s_inf[V_RESOLUTION];
     static double TT04LUT_EpiM_tau_s[V_RESOLUTION];
     static double TT04LUT_Endo_s_inf[V_RESOLUTION];
     static double TT04LUT_Endo_tau_s[V_RESOLUTION];
     static double TT04LUT_r_inf[V_RESOLUTION];
     static double TT04LUT_tau_r[V_RESOLUTION];     

           
     // lookup tables
     static double TT04LUT_rec_ipK[V_RESOLUTION];
     static double TT04LUT_NaCa_P1[V_RESOLUTION];
     static double TT04LUT_NaCa_P2[V_RESOLUTION];
     static double TT04LUT_NaK_P1[V_RESOLUTION];

     
     IBM_TT04_LUT();
     ~IBM_TT04_LUT();


     static void TT04LUT_Init();
     static void TT04LUT_InitParameters();
     static void TT04LUT_InitTable();

     static double Get_m_inf(int Vi){return(TT04LUT_m_inf[Vi]);};
     static double Get_h_inf(int Vi){return(TT04LUT_h_inf[Vi]);};
     static double Get_j_inf(int Vi){return(TT04LUT_j_inf[Vi]);};
     static double Get_Xr1_inf(int Vi){return(TT04LUT_Xr1_inf[Vi]);};
     static double Get_Xr2_inf(int Vi){return(TT04LUT_Xr2_inf[Vi]);};
     static double Get_Xs_inf(int Vi){return(TT04LUT_Xs_inf[Vi]);};
     static double Get_d_inf(int Vi){return(TT04LUT_d_inf[Vi]);};
     static double Get_f_inf(int Vi){return(TT04LUT_f_inf[Vi]);};
     static double Get_tau_m(int Vi){return(TT04LUT_tau_m[Vi]);};
     static double Get_tau_h(int Vi){return(TT04LUT_tau_h[Vi]);};
     static double Get_tau_j(int Vi){return(TT04LUT_tau_j[Vi]);};
     static double Get_tau_Xr1(int Vi){return(TT04LUT_tau_Xr1[Vi]);};
     static double Get_tau_Xr2(int Vi){return(TT04LUT_tau_Xr2[Vi]);};
     static double Get_tau_Xs(int Vi){return(TT04LUT_tau_Xs[Vi]);};
     static double Get_tau_d(int Vi){return(TT04LUT_tau_d[Vi]);};
     static double Get_tau_f(int Vi){return(TT04LUT_tau_f[Vi]);};
     static double Get_NaCa_P1(int Vi){return(TT04LUT_NaCa_P1[Vi]);};
     static double Get_NaCa_P2(int Vi){return(TT04LUT_NaCa_P2[Vi]);};
     static double Get_NaK_P1(int Vi){return(TT04LUT_NaK_P1[Vi]);};
     static double Get_CaL_P1(int Vi){return(TT04LUT_CaL_P1[Vi]);};
     static double Get_CaL_P2(int Vi){return(TT04LUT_CaL_P2[Vi]);};
     static double Get_r_inf(int Vi){return(TT04LUT_r_inf[Vi]);};
     static double Get_tau_r(int Vi){return(TT04LUT_tau_r[Vi]);};
     static double Get_EpiM_s_inf(int Vi){return(TT04LUT_EpiM_s_inf[Vi]);};
     static double Get_EpiM_tau_s(int Vi){return(TT04LUT_EpiM_tau_s[Vi]);};
     static double Get_Endo_s_inf(int Vi){return(TT04LUT_Endo_s_inf[Vi]);};
     static double Get_Endo_tau_s(int Vi){return(TT04LUT_Endo_tau_s[Vi]);};
     static double Get_rec_ipK(int Vi){return(TT04LUT_rec_ipK[Vi]);};
                                                                                                                    
};
                                                                                                                    
#endif
