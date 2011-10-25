/*
*/
#ifndef PARAMETER_INCLUDED
#include "IBM_tenTusscher04ParameterDefNLoad.hh"
#endif

class IBM_tenTusscher04_epiLUT
{
  public:
  
     // IBM parameter struct
     static char TT04LUT_EVfilename[255];// = "/bgpusr1/mr198/IBT/src/CellModel/ev_forceParameter/Flaim0910MutA/288_100.ev";
     static TTParameterStruct TT04LUT_ttps; //fps stands for TTParameterStruct
              
     // lookup tables
     static double TT04LUT_rec_ipK[V_RESOLUTION];
     static double TT04LUT_d_inf[V_RESOLUTION];
     static double TT04LUT_f_inf[V_RESOLUTION];
     static double TT04LUT_tau_m[V_RESOLUTION];
     static double TT04LUT_tau_h[V_RESOLUTION];
     static double TT04LUT_tau_j[V_RESOLUTION];
     static double TT04LUT_tau_Xr1[V_RESOLUTION];
     static double TT04LUT_tau_Xr2[V_RESOLUTION];
     static double TT04LUT_tau_Xs[V_RESOLUTION];
     static double TT04LUT_tau_r[V_RESOLUTION];
     static double TT04LUT_tau_s[V_RESOLUTION];
     static double TT04LUT_m_inf[V_RESOLUTION];
     static double TT04LUT_h_inf[V_RESOLUTION];
     static double TT04LUT_j_inf[V_RESOLUTION];
     static double TT04LUT_Xr1_inf[V_RESOLUTION];
     static double TT04LUT_Xr2_inf[V_RESOLUTION];
     static double TT04LUT_Xs_inf[V_RESOLUTION];
     static double TT04LUT_r_inf[V_RESOLUTION];
     static double TT04LUT_s_inf[V_RESOLUTION];
     static double TT04LUT_tau_d[V_RESOLUTION];
     static double TT04LUT_tau_f[V_RESOLUTION];
     static double TT04LUT_NaCa_P1[V_RESOLUTION];
     static double TT04LUT_NaCa_P2[V_RESOLUTION];
     static double TT04LUT_NaK_P1[V_RESOLUTION];
     static double TT04LUT_CaL_P1[V_RESOLUTION];
     static double TT04LUT_CaL_P2[V_RESOLUTION];
     
     IBM_tenTusscher04_epiLUT();
     ~IBM_tenTusscher04_epiLUT();


     static void TT04LUT_Init();
     static void TT04LUT_InitParameters();
     static void TT04LUT_InitTable();
            
     static double Get_rec_ipK(int Vi){return(TT04LUT_rec_ipK[Vi]);};
     static double Get_d_inf(int Vi){return(TT04LUT_d_inf[Vi]);};
     static double Get_f_inf(int Vi){return(TT04LUT_f_inf[Vi]);};
     static double Get_tau_m(int Vi){return(TT04LUT_tau_m[Vi]);};
     static double Get_tau_h(int Vi){return(TT04LUT_tau_h[Vi]);};
     static double Get_tau_j(int Vi){return(TT04LUT_tau_j[Vi]);};
     static double Get_tau_Xr1(int Vi){return(TT04LUT_tau_Xr1[Vi]);};
     static double Get_tau_Xr2(int Vi){return(TT04LUT_tau_Xr2[Vi]);};
     static double Get_tau_Xs(int Vi){return(TT04LUT_tau_Xs[Vi]);};
     static double Get_tau_r(int Vi){return(TT04LUT_tau_r[Vi]);};
     static double Get_tau_s(int Vi){return(TT04LUT_tau_s[Vi]);};
     static double Get_m_inf(int Vi){return(TT04LUT_m_inf[Vi]);};
     static double Get_h_inf(int Vi){return(TT04LUT_h_inf[Vi]);};
     static double Get_j_inf(int Vi){return(TT04LUT_j_inf[Vi]);};
     static double Get_Xr1_inf(int Vi){return(TT04LUT_Xr1_inf[Vi]);};
     static double Get_Xr2_inf(int Vi){return(TT04LUT_Xr2_inf[Vi]);};
     static double Get_Xs_inf(int Vi){return(TT04LUT_Xs_inf[Vi]);};
     static double Get_r_inf(int Vi){return(TT04LUT_r_inf[Vi]);};
     static double Get_s_inf(int Vi){return(TT04LUT_s_inf[Vi]);};
     static double Get_tau_d(int Vi){return(TT04LUT_tau_d[Vi]);};
     static double Get_tau_f(int Vi){return(TT04LUT_tau_f[Vi]);};
     static double Get_NaCa_P1(int Vi){return(TT04LUT_NaCa_P1[Vi]);};
     static double Get_NaCa_P2(int Vi){return(TT04LUT_NaCa_P2[Vi]);};
     static double Get_NaK_P1(int Vi){return(TT04LUT_NaK_P1[Vi]);};
     static double Get_CaL_P1(int Vi){return(TT04LUT_CaL_P1[Vi]);};
     static double Get_CaL_P2(int Vi){return(TT04LUT_CaL_P2[Vi]);};
                                                                                                                    
                                                                                                                    
};
                                                                                                                    
