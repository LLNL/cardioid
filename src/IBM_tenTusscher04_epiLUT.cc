/*
*/
#include "IBM_tenTusscher04_epiLUT.hh"



  char IBM_tenTusscher04_epiLUT::TT04LUT_EVfilename[255];
  TTParameterStruct IBM_tenTusscher04_epiLUT::TT04LUT_ttps; //fps stands for TTParameterStruct
     
  // lookup tables
  double IBM_tenTusscher04_epiLUT::TT04LUT_rec_ipK[V_RESOLUTION];
  double IBM_tenTusscher04_epiLUT::TT04LUT_d_inf[V_RESOLUTION];
  double IBM_tenTusscher04_epiLUT::TT04LUT_f_inf[V_RESOLUTION];
  double IBM_tenTusscher04_epiLUT::TT04LUT_tau_m[V_RESOLUTION];
  double IBM_tenTusscher04_epiLUT::TT04LUT_tau_h[V_RESOLUTION];
  double IBM_tenTusscher04_epiLUT::TT04LUT_tau_j[V_RESOLUTION];
  double IBM_tenTusscher04_epiLUT::TT04LUT_tau_Xr1[V_RESOLUTION];
  double IBM_tenTusscher04_epiLUT::TT04LUT_tau_Xr2[V_RESOLUTION];
  double IBM_tenTusscher04_epiLUT::TT04LUT_tau_Xs[V_RESOLUTION];
  double IBM_tenTusscher04_epiLUT::TT04LUT_tau_r[V_RESOLUTION];
  double IBM_tenTusscher04_epiLUT::TT04LUT_tau_s[V_RESOLUTION];
  double IBM_tenTusscher04_epiLUT::TT04LUT_m_inf[V_RESOLUTION];
  double IBM_tenTusscher04_epiLUT::TT04LUT_h_inf[V_RESOLUTION];
  double IBM_tenTusscher04_epiLUT::TT04LUT_j_inf[V_RESOLUTION];
  double IBM_tenTusscher04_epiLUT::TT04LUT_Xr1_inf[V_RESOLUTION];
  double IBM_tenTusscher04_epiLUT::TT04LUT_Xr2_inf[V_RESOLUTION];
  double IBM_tenTusscher04_epiLUT::TT04LUT_Xs_inf[V_RESOLUTION];
  double IBM_tenTusscher04_epiLUT::TT04LUT_r_inf[V_RESOLUTION];
  double IBM_tenTusscher04_epiLUT::TT04LUT_s_inf[V_RESOLUTION];
  double IBM_tenTusscher04_epiLUT::TT04LUT_tau_d[V_RESOLUTION];
  double IBM_tenTusscher04_epiLUT::TT04LUT_tau_f[V_RESOLUTION];
  double IBM_tenTusscher04_epiLUT::TT04LUT_NaCa_P1[V_RESOLUTION];
  double IBM_tenTusscher04_epiLUT::TT04LUT_NaCa_P2[V_RESOLUTION];
  double IBM_tenTusscher04_epiLUT::TT04LUT_NaK_P1[V_RESOLUTION];
  double IBM_tenTusscher04_epiLUT::TT04LUT_CaL_P1[V_RESOLUTION];
  double IBM_tenTusscher04_epiLUT::TT04LUT_CaL_P2[V_RESOLUTION];
                                                                                                                                            
                                                                                                                                                 
     
IBM_tenTusscher04_epiLUT::IBM_tenTusscher04_epiLUT() 
{

};


IBM_tenTusscher04_epiLUT::~IBM_tenTusscher04_epiLUT() 
{ 
};

void IBM_tenTusscher04_epiLUT::TT04LUT_Init()
{

//  printf("TT04LUT_Init() started\n");
//  fflush( stdout );

  TT04LUT_InitParameters();
  TT04LUT_InitTable();
    
//  printf("TT04LUT_Init() ended\n");
//  fflush( stdout );
        
    
}


void IBM_tenTusscher04_epiLUT::TT04LUT_InitParameters() 
{
// unused  int cellType = 0;

  TT04LUT_ttps.Rgas    = 8314.472;
  TT04LUT_ttps.Tx      = 310.0;
  TT04LUT_ttps.Faraday = 96485.3415;
  TT04LUT_ttps.K_o     = 5.4;
  TT04LUT_ttps.Ca_o    = 2.0;
  TT04LUT_ttps.Na_o    = 140.0;
  TT04LUT_ttps.Vc      = 0.016404;
  TT04LUT_ttps.Vsr     = 0.001094;
  TT04LUT_ttps.Bufc    = 0.15;
  TT04LUT_ttps.Kbufc   = 0.001;
  TT04LUT_ttps.Bufsr   = 10.0;
  TT04LUT_ttps.Kbufsr  = 0.3;
  TT04LUT_ttps.taufca  = 2.0;
  TT04LUT_ttps.taug    = 2.0;
  TT04LUT_ttps.Vmaxup  = 0.000425;
  TT04LUT_ttps.Kup     = 0.00025;
  TT04LUT_ttps.C       = 0.185;
  TT04LUT_ttps.g_Kr    = 0.096;
  TT04LUT_ttps.pKNa    = 0.03;
//  TT04LUT_ttps.g_Ks    = 0.245;
  TT04LUT_ttps.g_K1    = 5.405;
//  TT04LUT_ttps.g_to    = 0.073;
  TT04LUT_ttps.g_Na    = 14.838;
  TT04LUT_ttps.g_bNa   = 0.00029;
  TT04LUT_ttps.KmK     = 1.0;
  TT04LUT_ttps.KmNa    = 40.0;
  TT04LUT_ttps.knak    = 1.362;
  TT04LUT_ttps.g_CaL   = 0.000175;
  TT04LUT_ttps.g_bCa   = 0.000592;
  TT04LUT_ttps.kNaCa   = 1000;
  TT04LUT_ttps.KmNai   = 87.5;
  TT04LUT_ttps.KmCa    = 1.38;
  TT04LUT_ttps.ksat    = 0.1;
  TT04LUT_ttps.n       = 0.35;
  TT04LUT_ttps.g_pCa   = 0.825;
  TT04LUT_ttps.KpCa    = 0.0005;
  TT04LUT_ttps.g_pK    = 0.0146;
//  TT04LUT_ttps.s_inf_vHalf  = 28;
//  TT04LUT_ttps.tau_s_f1     = 1000;
//  TT04LUT_ttps.tau_s_slope1 = 1000;
//  TT04LUT_ttps.tau_s_vHalf1 = 67;
//  TT04LUT_ttps.tau_s_f2     = 0;
//  TT04LUT_ttps.tau_s_f3     = 8;

  // epi
       TT04LUT_ttps.g_to = 0.294;
       TT04LUT_ttps.s_inf_vHalf = 20;
       TT04LUT_ttps.tau_s_f1 = 85;
       TT04LUT_ttps.tau_s_slope1 = 320;
       TT04LUT_ttps.tau_s_vHalf1 = 45;
       TT04LUT_ttps.tau_s_f2 = 5;
       TT04LUT_ttps.tau_s_f3 = 3;
       TT04LUT_ttps.g_Ks = 0.245;


  // not in parameter file
  TT04LUT_ttps.m_Xr1_1 = -26;
  TT04LUT_ttps.m_Xr1_2 = 7;
  TT04LUT_ttps.a_Xr1_1 = -45;
  TT04LUT_ttps.a_Xr1_2 = 10;
  TT04LUT_ttps.b_Xr1_1 = -30;
  TT04LUT_ttps.b_Xr1_2 = 11.5;
  TT04LUT_ttps.K_Q10Xr1 = 1;
  
  TT04LUT_ttps.m_Xr2_1 = -88;
  TT04LUT_ttps.m_Xr2_2 = 24;
  TT04LUT_ttps.a_Xr2_1 = -60;
  TT04LUT_ttps.a_Xr2_2 = 20;
  TT04LUT_ttps.b_Xr2_1 = 60;
  TT04LUT_ttps.b_Xr2_2 = 20;
  TT04LUT_ttps.K_Q10Xr2 = 1;  

  //calculate remaining parameters
  TT04LUT_ttps.RToverF = TT04LUT_ttps.Rgas *( TT04LUT_ttps.Tx )/( TT04LUT_ttps.Faraday );
  TT04LUT_ttps.inverseRToverF =1/( TT04LUT_ttps.RToverF );
  TT04LUT_ttps.inverseVcF2C =(1/(2*( TT04LUT_ttps.Vc )*( TT04LUT_ttps.Faraday )))*( TT04LUT_ttps.C );
  TT04LUT_ttps.inverseVcFC =(1./(( TT04LUT_ttps.Vc )*( TT04LUT_ttps.Faraday )))*( TT04LUT_ttps.C );
  TT04LUT_ttps.VcdVsr = TT04LUT_ttps.Vc /( TT04LUT_ttps.Vsr );
  TT04LUT_ttps.Kupsquare = TT04LUT_ttps.Kup *( TT04LUT_ttps.Kup );
  TT04LUT_ttps.BufcPKbufc = TT04LUT_ttps.Bufc +( TT04LUT_ttps.Kbufc );
  TT04LUT_ttps.Kbufcsquare = TT04LUT_ttps.Kbufc *( TT04LUT_ttps.Kbufc );
  TT04LUT_ttps.Kbufc2 =2*( TT04LUT_ttps.Kbufc );
  TT04LUT_ttps.BufsrPKbufsr = TT04LUT_ttps.Bufsr +( TT04LUT_ttps.Kbufsr );
  TT04LUT_ttps.Kbufsrsquare = TT04LUT_ttps.Kbufsr *( TT04LUT_ttps.Kbufsr );
  TT04LUT_ttps.Kbufsr2 =2*( TT04LUT_ttps.Kbufsr );
  TT04LUT_ttps.KopKNaNao = TT04LUT_ttps.K_o +( TT04LUT_ttps.pKNa )*( TT04LUT_ttps.Na_o );
  TT04LUT_ttps.KmNai3 = TT04LUT_ttps.KmNai *( TT04LUT_ttps.KmNai )*( TT04LUT_ttps.KmNai );
  TT04LUT_ttps.Nao3 = TT04LUT_ttps.Na_o *( TT04LUT_ttps.Na_o )*( TT04LUT_ttps.Na_o );
}
                              
                              
                              

void IBM_tenTusscher04_epiLUT::TT04LUT_InitTable()
{

  for(double V=-RangeTabhalf+.0001; V<RangeTabhalf; V+=dDivisionTab) 
     {
       //V in mV
       const int Vi=(int)(DivisionTab*(RangeTabhalf+V)+.5);
                   
       
       const double rec_iNaK=(1./(1.+0.1245*exp(-0.1*V*( TT04LUT_ttps.inverseRToverF ))+0.0353*exp(-V*( TT04LUT_ttps.inverseRToverF ))));
       TT04LUT_NaK_P1[Vi]= TT04LUT_ttps.knak *( TT04LUT_ttps.K_o /( TT04LUT_ttps.K_o +( TT04LUT_ttps.KmK )))*rec_iNaK;
       TT04LUT_rec_ipK[Vi]=1./(1.+exp((25-V)/5.98));
       const double a_m=1./(1.+exp((-60.-V)/5.));
       const double b_m=0.1/(1.+exp((V+35.)/5.))+0.10/(1.+exp((V-50.)/200.));
       TT04LUT_tau_m[Vi]=a_m*b_m;
       TT04LUT_m_inf[Vi]=1./((1.+exp((-56.86-V)/9.03))*(1.+exp((-56.86-V)/9.03)));
       
       if (V>=-40.)
          {
            const double AH_1=0.;
            const double BH_1=(0.77/(0.13*(1.+exp(-(V+10.66)/11.1))));
            TT04LUT_tau_h[Vi]= 1.0/(AH_1+BH_1);
          }
       else
          {
            const double AH_2=(0.057*exp(-(V+80.)/6.8));
            const double BH_2=(2.7*exp(0.079*V)+(3.1e5)*exp(0.3485*V));
            TT04LUT_tau_h[Vi]=1.0/(AH_2+BH_2);
          }
       TT04LUT_h_inf[Vi]=1./((1.+exp((V+71.55)/7.43))*(1.+exp((V+71.55)/7.43)));
       
       if (V>=-40.)
          {
            const double AJ_1=0.;
            const double BJ_1=(0.6*exp((0.057)*V)/(1.+exp(-0.1*(V+32.))));
            TT04LUT_tau_j[Vi]= 1.0/(AJ_1+BJ_1);
          }
       else
          {
            const double AJ_2=(((-2.5428e4)*exp(0.2444*V)-(6.948e-6)
                             * exp(-0.04391*V))*(V+37.78)
                             / (1.+exp(0.311*(V+79.23))));
            const double BJ_2=(0.02424*exp(-0.01052*V)/(1.+exp(-0.1378*(V+40.14))));
            TT04LUT_tau_j[Vi]= 1.0/(AJ_2+BJ_2);
          }
       TT04LUT_j_inf[Vi]=TT04LUT_h_inf[Vi];
       TT04LUT_Xr1_inf[Vi]=1./(1.0+exp(( TT04LUT_ttps.m_Xr1_1 -V)/( TT04LUT_ttps.m_Xr1_2 )));
       const double a_Xr1=450./(1.+exp(( TT04LUT_ttps.a_Xr1_1 -V)/( TT04LUT_ttps.a_Xr1_2 )));
       const double b_Xr1=6./(1.+exp((V-( TT04LUT_ttps.b_Xr1_1 ))/( TT04LUT_ttps.b_Xr1_2 )));
       TT04LUT_tau_Xr1[Vi]=a_Xr1*b_Xr1*fabs( TT04LUT_ttps.K_Q10Xr1 );
       TT04LUT_Xr2_inf[Vi]=1./(1.+exp((V-( TT04LUT_ttps.m_Xr2_1 ))/( TT04LUT_ttps.m_Xr2_2 )));
       const double a_Xr2=3./(1.+exp(( TT04LUT_ttps.a_Xr2_1 -V)/( TT04LUT_ttps.a_Xr2_2 )));
       const double b_Xr2=1.12/(1.+exp((V-( TT04LUT_ttps.b_Xr2_1 ))/( TT04LUT_ttps.b_Xr2_2 )));
       TT04LUT_tau_Xr2[Vi]=a_Xr2*b_Xr2*fabs( TT04LUT_ttps.K_Q10Xr2 );
       TT04LUT_Xs_inf[Vi]=1./(1.+exp((-5.-V)/14.));
       const double a_Xs=1100./(sqrt(1.+exp((-10.-V)/6)));
       const double b_Xs=1./(1.+exp((V-60.)/20.));
       TT04LUT_tau_Xs[Vi]=a_Xs*b_Xs;
       TT04LUT_r_inf[Vi]=1./(1.+exp((20-V)/6.));
       TT04LUT_s_inf[Vi]=1./(1.+exp((V+( TT04LUT_ttps.s_inf_vHalf ))/5.));
       TT04LUT_tau_r[Vi]=9.5*exp(-(V+40.)*(V+40.)/1800.)+0.8;
       TT04LUT_tau_s[Vi]=( TT04LUT_ttps.tau_s_f1 )*exp(-(V+( TT04LUT_ttps.tau_s_vHalf1 ))*(V+( TT04LUT_ttps.tau_s_vHalf1 ))/( TT04LUT_ttps.tau_s_slope1 ))
                +( TT04LUT_ttps.tau_s_f2 )                                                         /(1.+exp((V-20.)/5.))
                +( TT04LUT_ttps.tau_s_f3 );
       
                                                
       
       TT04LUT_d_inf[Vi]=1./(1.+exp((-5-V)/7.5));
       const double a_d=1.4/(1.+exp((-35-V)/13))+0.25;
       const double b_d=1.4/(1.+exp((V+5)/5));
       const double c_d=1./(1.+exp((50-V)/20));
       TT04LUT_tau_d[Vi] =a_d*b_d+c_d;
       TT04LUT_f_inf[Vi]=1./(1.+exp((V+20)/7));
       TT04LUT_tau_f[Vi]=1125*exp(-(V+27)*(V+27)/240)+80+165/(1.+exp((25-V)/10));
       const double NaCaP1=( TT04LUT_ttps.kNaCa )*(1./(( TT04LUT_ttps.KmNai3 )+( TT04LUT_ttps.Nao3 )))
                          *(1./(( TT04LUT_ttps.KmCa )+( TT04LUT_ttps.Ca_o )))
                          *(1./(1.+( TT04LUT_ttps.ksat )*exp((( TT04LUT_ttps.n )-1.)
                          *V*( TT04LUT_ttps.inverseRToverF ))));
                                                                                                                               
       
       const double NaCaP2=exp(( TT04LUT_ttps.n )*V*( TT04LUT_ttps.inverseRToverF ))*( TT04LUT_ttps.Ca_o );
       const double NaCaP3=exp((( TT04LUT_ttps.n )-1)*V*( TT04LUT_ttps.inverseRToverF ))*( TT04LUT_ttps.Nao3 )*2.5;
       //INaCa=NaCaP1*(NaCaP2*Na_I^3-NaCaP3 TT04LUT_ttps.C a_i)=NaCa_P1[Vi]*Na_i^3-NaCa_P2[Vi] TT04LUT_ttps.C a_i , NaCa_P1[Vi]=NaCaP1*NaCaP2, NaCa_P2[
       TT04LUT_NaCa_P1[Vi]=NaCaP1*NaCaP2;
       TT04LUT_NaCa_P2[Vi]=NaCaP1*NaCaP3;
       const double CaLP1=2*V*( TT04LUT_ttps.inverseRToverF );
       const double CaLP2=2*( TT04LUT_ttps.g_CaL )*CaLP1*( TT04LUT_ttps.Faraday )/(exp(CaLP1)-1.);
       TT04LUT_CaL_P1[Vi]=CaLP2*exp(CaLP1);
       TT04LUT_CaL_P2[Vi]=CaLP2*-0.341*( TT04LUT_ttps.Ca_o );
       
       
     }// end for

}

