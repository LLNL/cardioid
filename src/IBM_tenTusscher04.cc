/*	File: IBM_tenTusscher04.cpp
    automatically created by ExtractParameterClass.pl - done by dw (20.04.2007)
    Institute of Biomedical Engineering, Universität Karlsruhe (TH)
    send comments to dw@ibt.uka.de	*/

#include "IBM_tenTusscher04.hh"

IBM_tenTusscher04::IBM_tenTusscher04(char *inputEVFilename, int cellType) 
{

  y_TT = (double*)malloc( sizeof(double) * TT04_STATE_VARIABLES);
  
#ifdef EULER
  dy_TT = (double*)malloc( sizeof(double) * TT04_STATE_VARIABLES);
#endif // EULER
#ifdef RK4
  yk = (double*)malloc( sizeof(double) * TT04_STATE_VARIABLES);
  k1 = (double*)malloc( sizeof(double) * TT04_STATE_VARIABLES);
  k2 = (double*)malloc( sizeof(double) * TT04_STATE_VARIABLES);
  k3 = (double*)malloc( sizeof(double) * TT04_STATE_VARIABLES);
  k4 = (double*)malloc( sizeof(double) * TT04_STATE_VARIABLES);
#endif // RK4

  if (
     (y_TT == NULL)
#ifdef EULER
   ||(dy_TT == NULL)
#endif  //EULER
#ifdef RK4
   ||(yk == NULL)||(k1 == NULL)||(k2 == NULL)||(k3 == NULL)||(k4 == NULL)
#endif // RK4
     )

     {
#ifdef USE_LOGLINES
       BegLogLine(1)
            << "Memory allocation failed for TT04_STATE_VARIABLES"
            << EndLogLine;
            
       fflush(stdout);
#endif
       exit(-1);
     }
                                                       
              

  Init(cellType);
//  printf("After INIT: ttps.RToverF %lf; ttps.K_o %lf; y_TT[tt_K_i] %lf\n", ttps.RToverF, ttps.K_o, y_TT[tt_K_i]);
//  fflush( stdout );
   

/*  
  strcpy(this->EVfilename, inputEVFilename);
  //cout << "setup evfile " << inputEVFilename << " " << this->EVfilename << endl;
  //cout << "init evfile " << this->EVfilename;
  loadTTParameter(ttps, y_TT, EVfilename);
  //cout << " ... done" << endl;
*/
  parameterCalculate();
                                              
//  printf("TT04 init done\n");
//  fflush( stdout );


};

IBM_tenTusscher04::~IBM_tenTusscher04() 
{ 
};


inline unsigned char IBM_tenTusscher04::getSpeed( double adVm) 
{
  return (unsigned char) 5;
};

void IBM_tenTusscher04::Init(int cellType) 
{
  my_cellType = cellType;

  y_TT[tt_V]    = -0.0862;
  y_TT[tt_m]    = 0.0; //(v(VT_m_init));
  y_TT[tt_h]    = 0.75;//(v(VT_h_init));
  y_TT[tt_j]    = 0.75;//(v(VT_j_init));
  y_TT[tt_xr1]  = 0;   //(v(VT_xr1_init));
  y_TT[tt_xr2]  = 1;   //(v(VT_xr2_init));
  y_TT[tt_xs]   = 0;   //(v(VT_xs_init));
  y_TT[tt_r]    = 0;   //(v(VT_r_init));
  y_TT[tt_s]    = 1;   //(v(VT_s_init));
  y_TT[tt_d]    = 0;   //(v(VT_d_init));
  y_TT[tt_f]    = 1;   //(v(VT_f_init));
  y_TT[tt_fCa]  = 1;   //(v(VT_fCa_init));
  y_TT[tt_g]    = 1;   //(v(VT_g_init));
  y_TT[tt_Ca_i] = 0.0002; //(v(VT_Cai_init));
  y_TT[tt_CaSR] = 0.2; //(v(VT_CaSR_init));
  y_TT[tt_Na_i] = 11.6;//(v(VT_Nai_init));
  y_TT[tt_K_i]  =fabs(138.3);

  ttps.Rgas    = 8314.472;
  ttps.Tx      = 310.0;
  ttps.Faraday = 96485.3415;
  ttps.K_o     = 5.4;
  ttps.Ca_o    = 2.0;
  ttps.Na_o    = 140.0;
  ttps.Vc      = 0.016404;
  ttps.Vsr     = 0.001094;
  ttps.Bufc    = 0.15;
  ttps.Kbufc   = 0.001;
  ttps.Bufsr   = 10.0;
  ttps.Kbufsr  = 0.3;
  ttps.taufca  = 2.0;
  ttps.taug    = 2.0;
  ttps.Vmaxup  = 0.000425;
  ttps.Kup     = 0.00025;
  ttps.C       = 0.185;
  ttps.g_Kr    = 0.096;
  ttps.pKNa    = 0.03;
//  ttps.g_Ks    = 0.245;
  ttps.g_K1    = 5.405;
//  ttps.g_to    = 0.073;
  ttps.g_Na    = 14.838;
  ttps.g_bNa   = 0.00029;
  ttps.KmK     = 1.0;
  ttps.KmNa    = 40.0;
  ttps.knak    = 1.362;
  ttps.g_CaL   = 0.000175;
  ttps.g_bCa   = 0.000592;
  ttps.kNaCa   = 1000;
  ttps.KmNai   = 87.5;
  ttps.KmCa    = 1.38;
  ttps.ksat    = 0.1;
  ttps.n       = 0.35;
  ttps.g_pCa   = 0.825;
  ttps.KpCa    = 0.0005;
  ttps.g_pK    = 0.0146;
//  ttps.s_inf_vHalf  = 28;
//  ttps.tau_s_f1     = 1000;
//  ttps.tau_s_slope1 = 1000;
//  ttps.tau_s_vHalf1 = 67;
//  ttps.tau_s_f2     = 0;
//  ttps.tau_s_f3     = 8;

  if (cellType == 0)  // endo
     {
       ttps.g_to = 0.073;
       ttps.s_inf_vHalf = 28;
       ttps.tau_s_f1 = 1000;
       ttps.tau_s_slope1 = 1000;
       ttps.tau_s_vHalf1 = 67;
       ttps.tau_s_f2 = 0;
       ttps.tau_s_f3 = 8;
       ttps.g_Ks = 0.245;
     }
  else if (cellType == 1) // mid
     {
       ttps.g_to = 0.294;
       ttps.s_inf_vHalf = 20;
       ttps.tau_s_f1 = 85;
       ttps.tau_s_slope1 = 320;
       ttps.tau_s_vHalf1 = 45;
       ttps.tau_s_f2 = 5;
       ttps.tau_s_f3 = 3;
       ttps.g_Ks = 0.062;
     }
  else // epi
     {
       ttps.g_to = 0.294;
       ttps.s_inf_vHalf = 20;
       ttps.tau_s_f1 = 85;
       ttps.tau_s_slope1 = 320;
       ttps.tau_s_vHalf1 = 45;
       ttps.tau_s_f2 = 5;
       ttps.tau_s_f3 = 3;
       ttps.g_Ks = 0.245;
     
     }


  // not in parameter file
  ttps.m_Xr1_1 = -26;
  ttps.m_Xr1_2 = 7;
  ttps.a_Xr1_1 = -45;
  ttps.a_Xr1_2 = 10;
  ttps.b_Xr1_1 = -30;
  ttps.b_Xr1_2 = 11.5;
  ttps.K_Q10Xr1 = 1;
  
  ttps.m_Xr2_1 = -88;
  ttps.m_Xr2_2 = 24;
  ttps.a_Xr2_1 = -60;
  ttps.a_Xr2_2 = 20;
  ttps.b_Xr2_1 = 60;
  ttps.b_Xr2_2 = 20;
  ttps.K_Q10Xr2 = 1;  




};



void IBM_tenTusscher04::parameterCalculate()
{

  ttps.RToverF = ttps.Rgas *( ttps.Tx )/( ttps.Faraday );
  ttps.inverseRToverF =1/( ttps.RToverF );
  ttps.inverseVcF2C =(1/(2*( ttps.Vc )*( ttps.Faraday )))*( ttps.C );
  ttps.inverseVcFC =(1./(( ttps.Vc )*( ttps.Faraday )))*( ttps.C );
  ttps.VcdVsr = ttps.Vc /( ttps.Vsr );
  ttps.Kupsquare = ttps.Kup *( ttps.Kup );
  ttps.BufcPKbufc = ttps.Bufc +( ttps.Kbufc );
  ttps.Kbufcsquare = ttps.Kbufc *( ttps.Kbufc );
  ttps.Kbufc2 =2*( ttps.Kbufc );
  ttps.BufsrPKbufsr = ttps.Bufsr +( ttps.Kbufsr );
  ttps.Kbufsrsquare = ttps.Kbufsr *( ttps.Kbufsr );
  ttps.Kbufsr2 =2*( ttps.Kbufsr );
  ttps.KopKNaNao = ttps.K_o +( ttps.pKNa )*( ttps.Na_o );
  ttps.KmNai3 = ttps.KmNai *( ttps.KmNai )*( ttps.KmNai );
  ttps.Nao3 = ttps.Na_o *( ttps.Na_o )*( ttps.Na_o );

};
  
                                                                                                                                        
double IBM_tenTusscher04::Calc(double dt, double V, double i_external)
{


  double svolt = V; //membrane voltage in mV
  double HT = dt; //timestep in ms
  const int Vi=(int)(DivisionTab*(RangeTabhalf+svolt)+.5); //array position


  // get values from LUTS
  if (my_cellType == 0)  // endo
     {
       rec_ipK = IBM_tenTusscher04_endoLUT::Get_rec_ipK(Vi);
       d_inf = IBM_tenTusscher04_endoLUT::Get_d_inf(Vi);
       f_inf = IBM_tenTusscher04_endoLUT::Get_f_inf(Vi);
       tau_m = IBM_tenTusscher04_endoLUT::Get_tau_m(Vi);
       tau_h = IBM_tenTusscher04_endoLUT::Get_tau_h(Vi);
       tau_j = IBM_tenTusscher04_endoLUT::Get_tau_j(Vi);
       tau_Xr1 = IBM_tenTusscher04_endoLUT::Get_tau_Xr1(Vi);
       tau_Xr2 = IBM_tenTusscher04_endoLUT::Get_tau_Xr2(Vi);
       tau_Xs = IBM_tenTusscher04_endoLUT::Get_tau_Xs(Vi);
       tau_r = IBM_tenTusscher04_endoLUT::Get_tau_r(Vi);
       tau_s = IBM_tenTusscher04_endoLUT::Get_tau_s(Vi);
       m_inf = IBM_tenTusscher04_endoLUT::Get_m_inf(Vi);
       h_inf = IBM_tenTusscher04_endoLUT::Get_h_inf(Vi);
       j_inf = IBM_tenTusscher04_endoLUT::Get_j_inf(Vi);
       Xr1_inf = IBM_tenTusscher04_endoLUT::Get_Xr1_inf(Vi);
       Xr2_inf = IBM_tenTusscher04_endoLUT::Get_Xr2_inf(Vi);
       Xs_inf = IBM_tenTusscher04_endoLUT::Get_Xs_inf(Vi);
       r_inf = IBM_tenTusscher04_endoLUT::Get_r_inf(Vi);
       s_inf = IBM_tenTusscher04_endoLUT::Get_s_inf(Vi);
       tau_d = IBM_tenTusscher04_endoLUT::Get_tau_d(Vi);
       tau_f = IBM_tenTusscher04_endoLUT::Get_tau_f(Vi);
       NaCa_P1 = IBM_tenTusscher04_endoLUT::Get_NaCa_P1(Vi);
       NaCa_P2 = IBM_tenTusscher04_endoLUT::Get_NaCa_P2(Vi);
       NaK_P1 = IBM_tenTusscher04_endoLUT::Get_NaK_P1(Vi);
       CaL_P1 = IBM_tenTusscher04_endoLUT::Get_CaL_P1(Vi);
       CaL_P2 = IBM_tenTusscher04_endoLUT::Get_CaL_P2(Vi);
     }
  else if (my_cellType == 1) // mid
     {
       rec_ipK = IBM_tenTusscher04_midLUT::Get_rec_ipK(Vi);
       d_inf = IBM_tenTusscher04_midLUT::Get_d_inf(Vi);
       f_inf = IBM_tenTusscher04_midLUT::Get_f_inf(Vi);
       tau_m = IBM_tenTusscher04_midLUT::Get_tau_m(Vi);
       tau_h = IBM_tenTusscher04_midLUT::Get_tau_h(Vi);
       tau_j = IBM_tenTusscher04_midLUT::Get_tau_j(Vi);
       tau_Xr1 = IBM_tenTusscher04_midLUT::Get_tau_Xr1(Vi);
       tau_Xr2 = IBM_tenTusscher04_midLUT::Get_tau_Xr2(Vi);
       tau_Xs = IBM_tenTusscher04_midLUT::Get_tau_Xs(Vi);
       tau_r = IBM_tenTusscher04_midLUT::Get_tau_r(Vi);
       tau_s = IBM_tenTusscher04_midLUT::Get_tau_s(Vi);
       m_inf = IBM_tenTusscher04_midLUT::Get_m_inf(Vi);
       h_inf = IBM_tenTusscher04_midLUT::Get_h_inf(Vi);
       j_inf = IBM_tenTusscher04_midLUT::Get_j_inf(Vi);
       Xr1_inf = IBM_tenTusscher04_midLUT::Get_Xr1_inf(Vi);
       Xr2_inf = IBM_tenTusscher04_midLUT::Get_Xr2_inf(Vi);
       Xs_inf = IBM_tenTusscher04_midLUT::Get_Xs_inf(Vi);
       r_inf = IBM_tenTusscher04_midLUT::Get_r_inf(Vi);
       s_inf = IBM_tenTusscher04_midLUT::Get_s_inf(Vi);
       tau_d = IBM_tenTusscher04_midLUT::Get_tau_d(Vi);
       tau_f = IBM_tenTusscher04_midLUT::Get_tau_f(Vi);
       NaCa_P1 = IBM_tenTusscher04_midLUT::Get_NaCa_P1(Vi);
       NaCa_P2 = IBM_tenTusscher04_midLUT::Get_NaCa_P2(Vi);
       NaK_P1 = IBM_tenTusscher04_midLUT::Get_NaK_P1(Vi);
       CaL_P1 = IBM_tenTusscher04_midLUT::Get_CaL_P1(Vi);
       CaL_P2 = IBM_tenTusscher04_midLUT::Get_CaL_P2(Vi);
     }
  else // epi
     {
       rec_ipK = IBM_tenTusscher04_epiLUT::Get_rec_ipK(Vi);
       d_inf = IBM_tenTusscher04_epiLUT::Get_d_inf(Vi);
       f_inf = IBM_tenTusscher04_epiLUT::Get_f_inf(Vi);
       tau_m = IBM_tenTusscher04_epiLUT::Get_tau_m(Vi);
       tau_h = IBM_tenTusscher04_epiLUT::Get_tau_h(Vi);
       tau_j = IBM_tenTusscher04_epiLUT::Get_tau_j(Vi);
       tau_Xr1 = IBM_tenTusscher04_epiLUT::Get_tau_Xr1(Vi);
       tau_Xr2 = IBM_tenTusscher04_epiLUT::Get_tau_Xr2(Vi);
       tau_Xs = IBM_tenTusscher04_epiLUT::Get_tau_Xs(Vi);
       tau_r = IBM_tenTusscher04_epiLUT::Get_tau_r(Vi);
       tau_s = IBM_tenTusscher04_epiLUT::Get_tau_s(Vi);
       m_inf = IBM_tenTusscher04_epiLUT::Get_m_inf(Vi);
       h_inf = IBM_tenTusscher04_epiLUT::Get_h_inf(Vi);
       j_inf = IBM_tenTusscher04_epiLUT::Get_j_inf(Vi);
       Xr1_inf = IBM_tenTusscher04_epiLUT::Get_Xr1_inf(Vi);
       Xr2_inf = IBM_tenTusscher04_epiLUT::Get_Xr2_inf(Vi);
       Xs_inf = IBM_tenTusscher04_epiLUT::Get_Xs_inf(Vi);
       r_inf = IBM_tenTusscher04_epiLUT::Get_r_inf(Vi);
       s_inf = IBM_tenTusscher04_epiLUT::Get_s_inf(Vi);
       tau_d = IBM_tenTusscher04_epiLUT::Get_tau_d(Vi);
       tau_f = IBM_tenTusscher04_epiLUT::Get_tau_f(Vi);
       NaCa_P1 = IBM_tenTusscher04_epiLUT::Get_NaCa_P1(Vi);
       NaCa_P2 = IBM_tenTusscher04_epiLUT::Get_NaCa_P2(Vi);
       NaK_P1 = IBM_tenTusscher04_epiLUT::Get_NaK_P1(Vi);
       CaL_P1 = IBM_tenTusscher04_epiLUT::Get_CaL_P1(Vi);
       CaL_P2 = IBM_tenTusscher04_epiLUT::Get_CaL_P2(Vi);
     }
                                                                                                                                                                                                            

//  printf("Vi %d; DivisionTab %d; RangeTabhalf %d; svolt %lf\n", Vi, DivisionTab, RangeTabhalf, svolt);
//  printf("tau_m %lf\n", tau_m);
//  fflush( stdout );
  
  //Needed to compute currents
  const double EK=( ttps.RToverF )*(log((( ttps.K_o )/y_TT[tt_K_i])));

//  printf("ttps.RToverF %lf; ttps.K_o %lf; y_TT[tt_K_i] %lf\n", ttps.RToverF, ttps.K_o, y_TT[tt_K_i]);
//  fflush( stdout );


  const double ENa=( ttps.RToverF )*(log((( ttps.Na_o )/y_TT[tt_Na_i])));
  const double EKs=( ttps.RToverF )*(log(( ttps.KopKNaNao )/(y_TT[tt_K_i]+( ttps.pKNa )*y_TT[tt_Na_i])));
  const double ECa=0.5*( ttps.RToverF )*(log((( ttps.Ca_o )/y_TT[tt_Ca_i])));
  const double AK1=0.1/(1.+exp(0.06*(svolt-EK-200)));
  const double BK1=(3.*exp(0.0002*(svolt-EK+100))+exp(0.1*(svolt-EK-10)))/(1.+exp(-0.5*(svolt-EK)));
  const double rec_iK1=AK1/(AK1+BK1);
  const double I_to=( ttps.g_to )*y_TT[tt_r]*y_TT[tt_s]*(svolt-EK);
  const double I_Kr=( ttps.g_Kr )*sqrt(( ttps.K_o )/5.4)*y_TT[tt_xr1]*y_TT[tt_xr2]*(svolt-EK);
  const double I_Ks=( ttps.g_Ks )*y_TT[tt_xs]*y_TT[tt_xs]*(svolt-EKs);
  const double I_K1=( ttps.g_K1 )*rec_iK1*(svolt-EK);
  const double I_pCa=( ttps.g_pCa )*y_TT[tt_Ca_i]/(( ttps.KpCa )+y_TT[tt_Ca_i]);
  const double I_pK=( ttps.g_pK )* rec_ipK  *(svolt-EK);
  const double I_bNa=( ttps.g_bNa )*(svolt-ENa);
  const double I_bCa=( ttps.g_bCa )*(svolt-ECa);

  //Compute scalar currents
  

  
  const double I_Na=( ttps.g_Na )*y_TT[tt_m]*y_TT[tt_m]*y_TT[tt_m]*y_TT[tt_h]*y_TT[tt_j]*(svolt-ENa);
  const double I_CaL=y_TT[tt_d]*y_TT[tt_f]*y_TT[tt_fCa]*( CaL_P1  *y_TT[tt_Ca_i]+ CaL_P2  );
  const double I_NaCa=( NaCa_P1  )*y_TT[tt_Na_i]*y_TT[tt_Na_i]*y_TT[tt_Na_i]-y_TT[tt_Ca_i]*( NaCa_P2  );
  const double I_NaK=( NaK_P1  )*(y_TT[tt_Na_i]/(y_TT[tt_Na_i]+( ttps.KmNa )));
  const double I_tot=I_Kr+I_Ks+I_K1+I_to+I_Na+I_bNa+I_CaL+I_bCa+I_NaK+I_NaCa+I_pCa+I_pK;

  //update concentrations
  const double Caisquare=y_TT[tt_Ca_i]*y_TT[tt_Ca_i];
  const double CaSRsquare=y_TT[tt_CaSR]*y_TT[tt_CaSR];
  const double CaCurrent=-(I_CaL+I_bCa+I_pCa-2*I_NaCa)*( ttps.inverseVcF2C );
  const double A=0.016464*CaSRsquare/(0.0625+CaSRsquare)+0.008232;
  const double I_rel=A*y_TT[tt_d]*y_TT[tt_g];
  const double I_leak=0.00008*(y_TT[tt_CaSR]-y_TT[tt_Ca_i]);
  const double SERCA=( ttps.Vmaxup )/(1.+(( ttps.Kupsquare )/Caisquare));
  const double HTCaSRCurrent=HT*(SERCA-I_rel-I_leak);
  const double CaCSQN=( ttps.Bufsr )*y_TT[tt_CaSR]/(y_TT[tt_CaSR]+( ttps.Kbufsr ));
  const double dCaSR=( ttps.VcdVsr )*HTCaSRCurrent;
  const double CaSum_1=CaCSQN+dCaSR+y_TT[tt_CaSR];
  const double bjsr=( ttps.BufsrPKbufsr )-CaSum_1;
  const double cjsr=( ttps.Kbufsr )*CaSum_1;
  y_TT[tt_CaSR]=(sqrt(bjsr*bjsr+4*cjsr)-bjsr)/2;
  const double CaBuf=( ttps.Bufc )*y_TT[tt_Ca_i]/(y_TT[tt_Ca_i]+( ttps.Kbufc ));
  const double dCai=HT*CaCurrent-HTCaSRCurrent;
  const double CaSum_2=CaBuf+dCai+y_TT[tt_Ca_i];
  const double bc=( ttps.BufcPKbufc )-CaSum_2;
  const double cc=( ttps.Kbufc )*CaSum_2;
  y_TT[tt_Ca_i]=(sqrt(bc*bc+4*cc)-bc)/2;
  const double dNai=-(I_Na+I_bNa+3*I_NaK+3*I_NaCa)*( ttps.inverseVcFC );
  y_TT[tt_Na_i]+=HT*dNai;
  const double dKi=-(I_to+I_Kr+I_Ks+I_K1-2*I_NaK+I_pK+i_external)*( ttps.inverseVcFC );
  y_TT[tt_K_i]+=HT*dKi;
  const double FCa_INF=(1./(1.+pow((y_TT[tt_Ca_i]*3076.923076923077),8))+0.1/(1.+exp(10000*y_TT[tt_Ca_i]-5))+0.20/(1.+exp(1250*y_TT[tt_Ca_i]-0.9375))+0.23)*.684931506849;
  const double G_INF=(y_TT[tt_Ca_i]<.00035?1./(1.+pow((y_TT[tt_Ca_i]*2857.142857142857),6)):1./(1.+pow((y_TT[tt_Ca_i]*2857.142857142857),16)));
  

  //update gates
  y_TT[tt_m] = m_inf  -(m_inf  -y_TT[tt_m])*exp(-HT/tau_m  );
  y_TT[tt_h] = h_inf  -(h_inf  -y_TT[tt_h])*exp(-HT/tau_h  );
  y_TT[tt_j] = j_inf  -(j_inf  -y_TT[tt_j])*exp(-HT/tau_j  );
  y_TT[tt_xr1] = Xr1_inf  -(Xr1_inf  -y_TT[tt_xr1])*exp(-HT/tau_Xr1  );
  y_TT[tt_xr2] = Xr2_inf  -(Xr2_inf  -y_TT[tt_xr2])*exp(-HT/tau_Xr2  );
  y_TT[tt_xs] = Xs_inf  -(Xs_inf  -y_TT[tt_xs])*exp(-HT/tau_Xs  );
  y_TT[tt_s]= s_inf  -(s_inf  -y_TT[tt_s])*exp(-HT/tau_s  );
  y_TT[tt_r]= r_inf  -(r_inf  -y_TT[tt_r])*exp(-HT/tau_r  );
  const double fcaold=y_TT[tt_fCa];
  const double gold=y_TT[tt_g];
  const double exptaufca=exp(-HT/( ttps.taufca ));
  const double exptaug=exp(-HT/( ttps.taug ));
  y_TT[tt_d] = d_inf  -(d_inf  -y_TT[tt_d])*exp(-HT/tau_d  ); 
  y_TT[tt_f] =f_inf  -(f_inf  -y_TT[tt_f])*exp(-HT/tau_f  ); 
  y_TT[tt_fCa] =FCa_INF-(FCa_INF-y_TT[tt_fCa])*exptaufca;
  y_TT[tt_g] =G_INF-(G_INF-y_TT[tt_g])*exptaug;
  if ((y_TT[tt_fCa]>fcaold) && ((svolt)>-60))
     {
       y_TT[tt_fCa]=fcaold;
     }  
       
  if ((y_TT[tt_g]>gold) && ((svolt)>-60))
     {
       y_TT[tt_g]=gold;
     }  
     
  const double dVmdt = (-I_tot);

  return(dVmdt); // should be in Volts/sec
}

void IBM_tenTusscher04::Print()//ostream &tempstr, double tArg,  double V) 
{  
/*
	tempstr<<tArg<<' '<<y_TT[tt_V]<<' '
	<<y_TT[tt_m]<<' '<<y_TT[tt_h]<<' '<<y_TT[tt_j]<<' '<<y_TT[tt_d]<<' '
	<<y_TT[tt_f]<<' '<<y_TT[tt_fCa]<<' '<<y_TT[tt_g]<<' '<<y_TT[tt_xr1]<<' '<<y_TT[tt_xr2]<<' '
	<<y_TT[tt_xs]<<' '<<y_TT[tt_r]<<' '<<y_TT[tt_s]<<' '<<y_TT[tt_Ca_i]<<' '<<y_TT[tt_CaSR]<<' '<<y_TT[tt_Na_i]<<' '<<y_TT[tt_K_i];
*/	
};

void IBM_tenTusscher04::LongPrint()//ostream &tempstr, double tArg,  double V) 
{
/*
    Print(tempstr, tArg, V);
      const double svolt=V*1000.0;
      constint Vi=(int)(DivisionTab*(RangeTabhalf+svolt)+.5);
	  const double EK=( ttps.RToverF)*(log((( ttps.K_o )/y_TT[tt_K_i])));
	  const double ENa=( ttps.RToverF )*(log((( ttps.Na_o )/y_TT[tt_Na_i])));
	  const double EKs=( ttps.RToverF )*(log(( ttps.KopKNaNao )/(y_TT[tt_K_i]+( ttps.pKNa )*y_TT[tt_Na_i])));
	  const double ECa=0.5*( ttps.RToverF )*(log((( ttps.Ca_o)/y_TT[tt_Ca_i])));
	  const double AK1=0.1/(1.+exp(0.06*(svolt-EK-200)));
	  const double BK1=(3.*exp(0.0002*(svolt-EK+100))+exp(0.1*(svolt-EK-10)))/(1.+exp(-0.5*(svolt-EK)));
	  const double rec_iK1=AK1/(AK1+BK1);
	  const double I_pCa=( ttps.g_pCa )*y_TT[tt_Ca_i]/(( ttps.KpCa )+y_TT[tt_Ca_i]);
	  const double I_pK=( ttps.g_pK )* rec_ipK  *(svolt-EK);
	  const double I_Na=( ttps.g_Na )*m*m*m*h*j*(svolt-ENa);
	  const double I_CaL=d*f*y_TT[tt_fCa]*( CaL_P1  *y_TT[tt_Ca_i]+ CaL_P2  );
	  const double I_NaCa=( NaCa_P1  )*y_TT[tt_Na_i]*y_TT[tt_Na_i]*y_TT[tt_Na_i]-y_TT[tt_Ca_i]*( NaCa_P2  );
	  const double I_NaK=( NaK_P1  )*(y_TT[tt_Na_i]/(y_TT[tt_Na_i]+( ttps.KmNa )));
	  const double I_to=( ttps.g_to )*y_TT[tt_r]*s*(svolt-EK);
	  const double I_Kr=( ttps.g_Kr )*sqrt(( ttps.K_o )/5.4)*y_TT[tt_xr1]*y_TT[tt_xr2]*(svolt-EK);
	  const double I_Ks=( ttps.g_Ks )*y_TT[tt_xs]*y_TT[tt_xs]*(svolt-EKs);
	  const double I_K1=( ttps.g_K1 )*rec_iK1*(svolt-EK);
	  const double I_bNa=( ttps.g_bNa )*(svolt-ENa);
	  const double I_bCa=( ttps.g_bCa )*(svolt-ECa);
	//  const double Caisquare=y_TT[tt_Ca_i]*y_TT[tt_Ca_i];
	  const double CaSRsquare=y_TT[tt_CaSR]*y_TT[tt_CaSR];
	//  const double CaCurrent=-(I_CaL+I_bCa+I_pCa-2*I_NaCa)*( ttps.inverseVcF2C );
	  const double A=0.016464*CaSRsquare/(0.0625+CaSRsquare)+0.008232;
	  const double I_rel=A*d*g;
	  const double I_leak=0.00008*(y_TT[tt_CaSR]-y_TT[tt_Ca_i]);
	  const double I_mem=I_Kr+I_Ks+I_K1+I_to+I_Na+I_bNa+I_CaL+I_bCa+I_NaK+I_NaCa+I_pCa+I_pK;
	tempstr<<' '<<I_Na/0.185<<' '<<I_CaL/0.185<<' '<<I_bCa/0.185<<' '<<I_pCa/0.185<<' '<<I_to/0.185<<' '<<I_Ks/0.185<<' '<<I_Kr/0.185<<' '<<I_K1/0.185
	    <<' '<<I_pK/0.185<<' '<<I_bNa/0.185<<' '<<I_NaK/0.185<<' '<<I_NaCa/0.185<<' '<<I_rel
	    <<' '<<I_leak <<' '<<I_mem << ' ';
*/
};
