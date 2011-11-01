/*	File: IBM_TT04.cc
*/

#include "IBM_TT04.hh"

IBM_TT04::IBM_TT04(char *inputEVFilename, int cellType) 
{

//    static bool tablesInitialized = false;
//    if ( !tablesInitialized )
//    {
//       tablesInitialized = true;
//       IBM_TT04_LUT::TT04LUT_Init();
//    }
   
   
  y_TT04 = (double*)malloc( sizeof(double) * TT04_STATE_VARIABLES);
  
#ifdef EULER
  dy_TT04 = (double*)malloc( sizeof(double) * TT04_STATE_VARIABLES);
#endif // EULER
#ifdef RK4
  yk = (double*)malloc( sizeof(double) * TT04_STATE_VARIABLES);
  k1 = (double*)malloc( sizeof(double) * TT04_STATE_VARIABLES);
  k2 = (double*)malloc( sizeof(double) * TT04_STATE_VARIABLES);
  k3 = (double*)malloc( sizeof(double) * TT04_STATE_VARIABLES);
  k4 = (double*)malloc( sizeof(double) * TT04_STATE_VARIABLES);
#endif // RK4

  if (
     (y_TT04 == NULL)
#ifdef EULER
   ||(dy_TT04 == NULL)
#endif  //EULER
#ifdef RK4
   ||(yk == NULL)||(k1 == NULL)||(k2 == NULL)||(k3 == NULL)||(k4 == NULL)
#endif // RK4
     )

     {
       printf("Memory allocation failed for state variables\n");
       exit(-1);
     }
              

  Init(cellType);

};

IBM_TT04::~IBM_TT04() 
{ 
};


void IBM_TT04::Init(int cellType) 
{
  my_cellType = cellType;
  
  // no given value in CellML - initialize to 0
  a_K1 = 0.0;
  b_K1 = 0.0;
  XK1_inf = 0.0;

  Xr1_inf = 0.0;
  tau_Xr1 = 0.0;
  Xr2_inf = 0.0;
  tau_Xr2 = 0.0;
  
  Xs_inf = 0.0;
  tau_Xs = 0.0;    
  
  m_inf = 0.0;
  tau_m = 0.0;
  h_inf = 0.0;
  tau_h = 0.0;
  j_inf = 0.0;
  tau_j = 0.0;

  d_inf = 0.0;
  tau_d = 0.0;
  f_inf = 0.0;
  tau_f = 0.0;
  alpha_fCa = 0.0;
  beta_fCa = 0.0;
  gamma_fCa = 0.0;
  fCa_inf = 0.0;
  tau_fCa = 2.0; // in [ms]
  tau_fCa_reciproc = 1 / tau_fCa; // in [ms]
  CaL_P1 = 0.0;
  CaL_P2 = 0.0;


  g_inf = 0.0;
  tau_g = 2.0; // in [ms]
  tau_g_reciproc = 1/tau_g;

  tau_r = 0.0;
  tau_s = 0.0;
  r_inf = 0.0;
  s_inf = 0.0;
  NaCa_P1 = 0.0;
  NaCa_P2 = 0.0;
  NaK_P1 = 0.0;

  rec_ipK = 0.0;

  Ca_ibufc = 0.0;
  Ca_srbufsr = 0.0;

  
  
  // initialize currents
  I_K1 = 0.0;
  I_to = 0.0;
  I_Kr = 0.0;
  I_Ks = 0.0;
  I_CaL = 0.0;
  I_NaK = 0.0;
  I_Na = 0.0;
  I_bNa = 0.0;
  I_NaCa = 0.0;
  I_bCa = 0.0;
  I_pK = 0.0;
  I_pCa = 0.0;
  I_leak = 0.0;
  I_up = 0.0;
  I_rel = 0.0;
  I_ax = 0.0;
  Iion = 0.0;
  
  ENa = 0.0;
  EK = 0.0;
  ECa = 0.0;
  EKs = 0.0;

  
/*
   Initial Conditions
*/   
  y_TT04[tt04_V]    = -86.2;  // -86.2 mV
  y_TT04[tt04_Nai] = 11.6;
  y_TT04[tt04_Ki]  = 138.3;
  y_TT04[tt04_Cai] = 0.0002; 
  y_TT04[tt04_xr1]  = 0.0;   
  y_TT04[tt04_xr2]  = 1.0;   
  y_TT04[tt04_xs]   = 0.0;   
  y_TT04[tt04_m]    = 0.0;
  y_TT04[tt04_h]    = 0.75;
  y_TT04[tt04_j]    = 0.75;
  y_TT04[tt04_d]    = 0.0;   
  y_TT04[tt04_f]    = 1.0;   
  y_TT04[tt04_fCa]  = 1.0;   
  y_TT04[tt04_s]    = 1.0;   
  y_TT04[tt04_r]    = 0.0;   
  y_TT04[tt04_CaSR] = 0.2; 
  y_TT04[tt04_g]    = 1.0;  
  
  
  // constants CellML tenTusscher et al. AJP Heart Circ Physiol 2004
  tt04ps.R = 8314.472;         // gas constant                     - J*K^(-1)*mol^(-1)   CellML Note: value different from TT04 publication
  tt04ps.T = 310.0;            // temperature                      - K
  tt04ps.F = 96485.3415;          // Faraday constant                 - C/mmol           CellML Note: value different from TT04 publication
  tt04ps.Cm = 0.185;               // Cell capacitance per unit surface area  - microF/cm^2
  tt04ps.VC = 0.016404;           // Cytoplasmic volume               - microm^3         CellML Note: value different from TT04 publication
  tt04ps.pKNa = 0.03;          // relative IKs permeability to Na+
  tt04ps.Nao = 140.0;            // Extracellular Na+ concentration  - mM
  tt04ps.Ko = 5.4;             // Extracellular K+ concentration   - mM
  tt04ps.Cao = 2.0;              // Extracellular Ca2+ concentration - mM
  tt04ps.GK1 = 5.405;          // Maximal IK1 conductance          - nS/pF
  tt04ps.GKr = 0.096;          // Maximal IKr conductance          - nS/pF
  tt04ps.GKs_epi_endo = 0.245; // Maximal IKs (epi/endo) conductance - nS/pF
  tt04ps.GKs_M = 0.062;        // Maximal IKs (M) conductance        - nS/pF
  tt04ps.GKs = 0.0;  // set later to cell type value
  tt04ps.GNa = 14.838;         // Maximal INa conductance          - nS/pF
  tt04ps.GbNa = 0.00029;       // Maximal IbNa conductance               - nS/pF
  tt04ps.GCaL = 0.000175;      // = 1.75*10^(-4) != {0.10662224 = 1.75^(-4)} Maximal ICaL conductance - cm^3*microF^(-1)*s^(-1)
  tt04ps.GbCa = 0.000592;      // Maximal IbCa conductance               - nS/pF
  tt04ps.Gto_epi_M = 0.294;    // Maximal Gto (epi/M) conductance  - nS/pF
  tt04ps.Gto_endo = 0.073;     // Maximal Gto (endo) conductance   - nS/pF
  tt04ps.Gto = 0.0;  // will be set later according to cell type
  tt04ps.PNaK = 1.362;         // Maximal INaK                           - pA/pF
  tt04ps.KmK = 1.0;              // Ko half-saturation constant of INaK    - mM
  tt04ps.KmNa = 40.0;            // Nai half-saturation constant of INaK   - mM
  tt04ps.kNaCa = 1000;         // Maximal INaCa                          - pA/pF
  tt04ps.ksat = 0.1;           // Saturation factor for INaCa
  tt04ps.alpha = 2.5;          // Factor enhancing outward nature of INa CA
  tt04ps.gamma = 0.35;         // Voltage dependence parameter of INaCa 
  tt04ps.KmCa = 1.38;          // Cai half saturation constant for INaCa - mM
  tt04ps.KmNai = 87.5;         // Nai half saturation constant for INaCa - mM
  tt04ps.GpCa = 0.825;         // Maximal IpCa conductance               - nS/pF    CellML Note: the TT04 publication has wrong number (0.025) cited
  tt04ps.KpCa = 0.0005;        // Cai half-saturation constant of IpCa   - mM
  tt04ps.GpK = 0.0146;         // Maximal IpK conductance                - nS/pF
  tt04ps.arel = 0.016464;        // Maximal CaSR-dependent Irel            - mM/s   CellML Note: value different from TT04 publication
  tt04ps.brel = 0.25;          // CaSR half-saturation constant of Irel  - mM
  tt04ps.crel = 0.008232;         // Maximal CaSR-independent Irel          - mM/s  CellML Note: value different from TT04 publication
  tt04ps.Kup = 0.00025;        // Half-saturation constant of Iup        - mM
  tt04ps.Vleak = 0.00008;      // Maximal Ileak                          - ms^(-1)
  tt04ps.Vmaxup = 0.000425;    // Maximal Iup                            - mM/ms
  tt04ps.Bufc = 0.15;          // Total cytoplasmic buffer concentration - mM
  tt04ps.Kbufc = 0.001;        // Cai half-saturation constant for cytoplasmic buffer concentraton - mM
  tt04ps.Bufsr = 10.0;           // Total sarcoplasmic buffer concentration - mM
  tt04ps.Kbufsr = 0.3;         // CaSR half-saturation constant for sarcoplasmic buffer - mM
  tt04ps.VSR = 0.001094;           // Sarcoplasmic reticulum volume    - microm^3  CellML Note: value different from TT04 publication
  
/*  
  //not used in single cell computation

  tt04ps.S = 0.2;              // surface to volume ratio          - microm^(-1)
  tt04ps.rho = 162;            // cellular resistivity             - Ohm*cm

*/  
  
  Kupsquare = tt04ps.Kup * tt04ps.Kup;
  brelsquare = tt04ps.brel * tt04ps.brel;

  // calculate parameters
  tt04ps.RToverF = (tt04ps.R * tt04ps.T )/( tt04ps.F );
  tt04ps.inverseRToverF = 1 / tt04ps.RToverF;
  tt04ps.KopKNaNao = tt04ps.Ko +( tt04ps.pKNa )*( tt04ps.Nao );
  tt04ps.CaL_4FF_RT = 4 * tt04ps.inverseRToverF * tt04ps.F;
  tt04ps.inverseVcF = (-1.0)*(1./(( tt04ps.VC )*( tt04ps.F)));
  tt04ps.inverse2VcF = 0.5 * tt04ps.inverseVcF;

  // only used in LUT computation but initialized nevertheless here as well  
  tt04ps.KmNai3 = tt04ps.KmNai *( tt04ps.KmNai )*( tt04ps.KmNai );
  tt04ps.Nao3   = tt04ps.Nao *( tt04ps.Nao )*( tt04ps.Nao );
  


  if (cellType == 0)  // endo
     {
       tt04ps.Gto = tt04ps.Gto_endo;
       tt04ps.GKs = tt04ps.GKs_epi_endo;
     }
  else if (cellType == 1) // mid
     {
       tt04ps.Gto = tt04ps.Gto_epi_M;
       tt04ps.GKs = tt04ps.GKs_M;
     }
  else if (cellType == 2)// epi
     {
       tt04ps.Gto = tt04ps.Gto_epi_M;
       tt04ps.GKs = tt04ps.GKs_epi_endo;
     }
  else
     {
       printf("Cell type out of scope\n");
       exit(-1);
     }

};

  
 
                                                                                                                                       
//double IBM_TT04::Calc(double dt,  double Vm,  double i_stim) // i_stim in CellML = 52 [pA/pF]
double IBM_TT04::Calc(double dt_ms,  double Vm_mV,  double i_stim) // i_stim in CellML = 52 [pA/pF]
{

//  double Vm_mV = Vm*1000; //membrane voltage in mV
//  double dt_ms = dt*1000; //timestep in ms
  const int Vi=(int)(DivisionTab*(RangeTabhalf+Vm_mV)+.5); //array position
                                                                                                                                                                                                            

//  printf("Vi %d; DivisionTab %d; RangeTabhalf %d; Vm_mV %lf\n", Vi, DivisionTab, RangeTabhalf, Vm_mV);
//  fflush( stdout );
  
  // Reversal Potentials
  // Eq. 25 - TT04
  ENa=( tt04ps.RToverF )*(log((( tt04ps.Nao )/y_TT04[tt04_Nai])));
  EK =( tt04ps.RToverF )*(log((( tt04ps.Ko  )/y_TT04[tt04_Ki] )));
  // Note: the factor 0.5 is not given in publication but in CellML
  ECa=(0.5 * tt04ps.RToverF)*(log((( tt04ps.Cao )/y_TT04[tt04_Cai])));
  // Eq. 26 - TT04
  EKs=( tt04ps.RToverF )*(log(( tt04ps.KopKNaNao )
                                     /(y_TT04[tt04_Ki]+( tt04ps.pKNa )*y_TT04[tt04_Nai])));
				     

  // Fast Na+ Current  
  // Eq. 4 = Eq. 27 - TT04
  I_Na=( tt04ps.GNa )*y_TT04[tt04_m]*y_TT04[tt04_m]*y_TT04[tt04_m]
                     *y_TT04[tt04_h]*y_TT04[tt04_j]*(Vm_mV - ENa);
  
  // Eq. 28 - 39 TT04
  m_inf = IBM_TT04_LUT::Get_m_inf(Vi);
  h_inf = IBM_TT04_LUT::Get_h_inf(Vi);
  j_inf = IBM_TT04_LUT::Get_j_inf(Vi);
  tau_m = IBM_TT04_LUT::Get_tau_m(Vi);
  tau_h = IBM_TT04_LUT::Get_tau_h(Vi);
  tau_j = IBM_TT04_LUT::Get_tau_j(Vi);


  // L-type Ca2+ current
  // Eq. 40 TT04
  I_CaL = y_TT04[tt04_d]*y_TT04[tt04_f]*y_TT04[tt04_fCa]
	* ((y_TT04[tt04_Cai] * CaL_P1) + CaL_P2);

  CaL_P1 = IBM_TT04_LUT::Get_CaL_P1(Vi);
  CaL_P2 = IBM_TT04_LUT::Get_CaL_P2(Vi);
  
  // Eq. 41- 47 - TT04
  d_inf = IBM_TT04_LUT::Get_d_inf(Vi);
  f_inf = IBM_TT04_LUT::Get_f_inf(Vi);
  tau_d = IBM_TT04_LUT::Get_tau_d(Vi);
  tau_f = IBM_TT04_LUT::Get_tau_f(Vi);
  
  // Eq. 48 - 51 - TT04
  // these are [Cai] dependent - hence not in a LUT with [V] base
  // Eq. 48 - TT04
  // Jim & Chanhoan: what is the most optimal way to compute ^8??
  // Maybe we need to think about [Cai] in a LUT similar to [V]
  // Note: f_Ca seemed to be calculated differently before...
  double a_fCa_denom1 = (y_TT04[tt04_Cai]/0.000325);
  double a_fCa_denom_pow8 = a_fCa_denom1 * a_fCa_denom1
                          * a_fCa_denom1 * a_fCa_denom1
			  * a_fCa_denom1 * a_fCa_denom1
                          * a_fCa_denom1 * a_fCa_denom1;
  
  alpha_fCa = 1/(1 + a_fCa_denom_pow8);
  // Eq. 49 - TT04
  beta_fCa = 0.1/(1 + exp((y_TT04[tt04_Cai]-0.0005)/0.0001));
  // Eq. 50 - TT04
  gamma_fCa = 0.2/(1 + exp((y_TT04[tt04_Cai]-0.00075)/0.0008));
  // Eq. 51 - TT04
  fCa_inf = (alpha_fCa + beta_fCa + gamma_fCa + 0.23)/1.46;
  // Eq. 52 - TT05
  // tau_fCa = 2.0; // constant - initialized above
  
  // Eq. 53 - 54 - TT04 see below where state variables are updated.


  // Transient outward current       
  // Eq. 55 - TT04
  I_to = tt04ps.Gto * y_TT04[tt04_r] * y_TT04[tt04_s] * (Vm_mV - EK);
  // Eq. 55 - TT04
  r_inf = IBM_TT04_LUT::Get_r_inf(Vi);
  // Eq. 56 - TT04
  tau_r = IBM_TT04_LUT::Get_tau_r(Vi);
  
  // Eq. 58 - 61 - TT04
  if ((my_cellType == 1)||(my_cellType == 2)) // mid or epi
     {
       s_inf = IBM_TT04_LUT::Get_EpiM_s_inf(Vi);
       tau_s = IBM_TT04_LUT::Get_EpiM_tau_s(Vi);
     }
  else if (my_cellType == 0) // endo
     {
       s_inf = IBM_TT04_LUT::Get_Endo_s_inf(Vi);
       tau_s = IBM_TT04_LUT::Get_Endo_tau_s(Vi);
     }
  else 
     {
       printf("Unknown cell type\n");
       exit(-1);
     }
  

  // Slow delayed rectifier current
  // Eq. 62 - TT04
  I_Ks = tt04ps.GKs * y_TT04[tt04_xs] * y_TT04[tt04_xs] * (Vm_mV - EKs);
  // Eq. 63 - 66 - TT04
  Xs_inf = IBM_TT04_LUT::Get_Xs_inf(Vi);
  tau_Xs = IBM_TT04_LUT::Get_tau_Xs(Vi);

  // Rapid delayed rectifier current
  // Eq. 67 - TT04
  I_Kr = tt04ps.GKr * sqrt((tt04ps.Ko/5.4)) * y_TT04[tt04_xr1] * y_TT04[tt04_xr2] * (Vm_mV - EK);  
  // Eq. 68 - 75 - TT04  
  Xr1_inf = IBM_TT04_LUT::Get_Xr1_inf(Vi);
  tau_Xr1 = IBM_TT04_LUT::Get_tau_Xr1(Vi);
  Xr2_inf = IBM_TT04_LUT::Get_Xr2_inf(Vi);
  tau_Xr2 = IBM_TT04_LUT::Get_tau_Xr2(Vi);


  // Inward rectifier K+ current
  // Eq. 76 - TT04
  // Note the sqrt could be stored in tmp and computed only once - see eq. 67
  I_K1 = tt04ps.GK1 * sqrt((tt04ps.Ko/5.4)) * XK1_inf * (Vm_mV - EK);
  // Eq. 77 - TT04
  a_K1=0.1/(1.+exp(0.06*(Vm_mV-EK-200)));
  b_K1=((3.0*exp(0.0002*(Vm_mV-EK+100))) + exp(0.1*(Vm_mV-EK-10)))
      /(1.+exp(-0.5*(Vm_mV-EK)));

  XK1_inf = a_K1/(a_K1+b_K1);

  // Na+/Ca2+ exchanger current
  // Eq.80 - TT04
  I_NaCa = (NaCa_P1 * y_TT04[tt04_Nai] * y_TT04[tt04_Nai] * y_TT04[tt04_Nai])
         + (y_TT04[tt04_Cai] * NaCa_P2);  // the minus is in NaCa_P2
  NaCa_P1 = IBM_TT04_LUT::Get_NaCa_P1(Vi);
  NaCa_P2 = IBM_TT04_LUT::Get_NaCa_P2(Vi);


  // Na+/K+ pump current
  // Eq. 81 - TT04
  I_NaK = (NaK_P1 * y_TT04[tt04_Nai])
        / (y_TT04[tt04_Nai] + tt04ps.KmNa);
  NaK_P1 = IBM_TT04_LUT::Get_NaK_P1(Vi);
  // Eq. 82 - TT04
  I_pCa = (tt04ps.GpCa * y_TT04[tt04_Cai])/( tt04ps.KpCa + y_TT04[tt04_Cai]);
  // Eq. 83 - TT04
  I_pK  = tt04ps.GpK  * rec_ipK  *(Vm_mV - EK);
  rec_ipK = IBM_TT04_LUT::Get_rec_ipK(Vi);


  // Background currents
  // Eq. 84 - TT04
  I_bNa = tt04ps.GbNa * (Vm_mV - ENa);  
  // Eq. 85 - TT04
  I_bCa = tt04ps.GbCa * (Vm_mV - ECa);

  // Calcium dynamics
  
  double Caisquare  = y_TT04[tt04_Cai] *y_TT04[tt04_Cai];
  double CaSRsquare = y_TT04[tt04_CaSR]*y_TT04[tt04_CaSR];

  
  // Eq. 86 - TT04
  I_leak = tt04ps.Vleak * (y_TT04[tt04_CaSR]-y_TT04[tt04_Cai]);
  // Eq. 87 - TT04
  I_up = tt04ps.Vmaxup / (1.0+(Kupsquare /Caisquare));
  // Eq. 87 - TT04
  I_rel = (((tt04ps.arel * CaSRsquare)/(brelsquare + CaSRsquare)) + tt04ps.crel)
        * y_TT04[tt04_d] * y_TT04[tt04_g];
  
  // Eq. 89 - TT04
  if (y_TT04[tt04_Cai] <= 0.00035)
     {
       double g_infDenom6 = (y_TT04[tt04_Cai]/0.00035);
       g_infDenom6 = g_infDenom6 * g_infDenom6 * g_infDenom6;
       g_infDenom6 *= g_infDenom6;
       g_inf = 1/ (1 + g_infDenom6);
     }
  else
     {
       double g_infDenom16 = (y_TT04[tt04_Cai]/0.00035);
       g_infDenom16 *= g_infDenom16;  //^2
       g_infDenom16 *= g_infDenom16;  //^4
       g_infDenom16 *= g_infDenom16;  //^16

       g_inf = 1/ (1 + g_infDenom16);
     }

  // Eq. 90 - TT04
  // tau_g = 2.0; see initialization
  // Eq. 91 - 92 - TT04 see below where state variables are updated.
  
  // Eq. 93 - TT04
  Ca_ibufc = 1 / (1 + ((tt04ps.Bufc * tt04ps.Kbufc)
                       /((y_TT04[tt04_Cai] + tt04ps.Kbufc)*(y_TT04[tt04_Cai] + tt04ps.Kbufc))));
  // Eq. 95 - TT04
  Ca_srbufsr = 1 / (1 + (tt04ps.Bufsr * tt04ps.Kbufsr)
                       /((y_TT04[tt04_CaSR] + tt04ps.Kbufsr)*(y_TT04[tt04_CaSR] + tt04ps.Kbufsr)));
  
  
  
  // UPDATE STATE VARIABLES 
  //update gates - exact solution would be:
  y_TT04[tt04_m] = (m_inf  -(m_inf  -y_TT04[tt04_m])*exp(-dt_ms/tau_m ));
  y_TT04[tt04_h] = h_inf  -(h_inf  -y_TT04[tt04_h])*exp(-dt_ms/tau_h  );
  y_TT04[tt04_j] = j_inf  -(j_inf  -y_TT04[tt04_j])*exp(-dt_ms/tau_j  );
  y_TT04[tt04_d] = d_inf  -(d_inf  -y_TT04[tt04_d])*exp(-dt_ms/tau_d  ); 
  y_TT04[tt04_f] =f_inf  -(f_inf  -y_TT04[tt04_f])*exp(-dt_ms/tau_f  ); 
  y_TT04[tt04_r]= r_inf  -(r_inf  -y_TT04[tt04_r])*exp(-dt_ms/tau_r  );
  y_TT04[tt04_s]= s_inf  -(s_inf  -y_TT04[tt04_s])*exp(-dt_ms/tau_s  );
  y_TT04[tt04_xs] = Xs_inf  -(s_inf  -y_TT04[tt04_xs])*exp(-dt_ms/tau_Xs  );
  y_TT04[tt04_xr1] = Xr1_inf  -(Xr1_inf  -y_TT04[tt04_xr1])*exp(-dt_ms/tau_Xr1  );
  y_TT04[tt04_xr2] = Xr2_inf  -(Xr2_inf  -y_TT04[tt04_xr2])*exp(-dt_ms/tau_Xr2  );

  //USING EXPLICIT EULER  
  // because we potentially truncate exp taylor series expantion by using LUT, we might as well just take Eulers
  // method
/*  
  y_TT04[tt04_m] += (dt_ms * ((m_inf  - y_TT04[tt04_m]) / tau_m ));
  y_TT04[tt04_h] += (dt_ms * ((h_inf  - y_TT04[tt04_h]) / tau_h ));
  y_TT04[tt04_j] += (dt_ms * ((j_inf  - y_TT04[tt04_j]) / tau_j ));
  y_TT04[tt04_d] += (dt_ms * ((d_inf  - y_TT04[tt04_d]) / tau_d )); 
  y_TT04[tt04_f] += (dt_ms * ((f_inf  - y_TT04[tt04_f]) / tau_f )); 
  y_TT04[tt04_r]   += (dt_ms * ((r_inf    - y_TT04[tt04_r])  /tau_r));
  y_TT04[tt04_s]   += (dt_ms * ((s_inf    - y_TT04[tt04_s])  /tau_s));
  y_TT04[tt04_xs]  += (dt_ms * ((Xs_inf   - y_TT04[tt04_xs]) /tau_Xs));
  y_TT04[tt04_xr1] += (dt_ms * ((Xr1_inf  - y_TT04[tt04_xr1])/tau_Xr1));
  y_TT04[tt04_xr2] += (dt_ms * ((Xr2_inf  - y_TT04[tt04_xr2])/tau_Xr2));
*/
  
  // Eq. 53 - 54 - TT04 see below where state variables are updated.
  if ((fCa_inf>y_TT04[tt04_fCa]) && ((Vm_mV)>(-60.0)))
     {
       y_TT04[tt04_fCa] = 0;  // k = 0
     }  
  else
     {
       y_TT04[tt04_fCa] += dt_ms * (fCa_inf - y_TT04[tt04_fCa]) * tau_fCa_reciproc; // k = 1
     }

  // Eq. 91 - 92 - TT04
  if ((g_inf > y_TT04[tt04_g]) && (Vm_mV > -60.0))
     {
       y_TT04[tt04_g] = 0;  // k = 0
     }
  else
     {
       y_TT04[tt04_g] += dt_ms * ((g_inf - y_TT04[tt04_g]) * tau_g_reciproc);  // k = 1
     }



  // update ion concentrations
  // Eq. 94 - TT04
  // the minus is in the tt04ps.inverseVcF factor
  // from paper: Caitotal = Caibufc + Cai
/*
  y_TT04[tt04_Cai] += dt_ms *  
                            (Ca_ibufc * (  (tt04ps.Cm * tt04ps.inverse2VcF * (I_CaL + I_bCa + I_pCa - (2 * I_NaCa)))
                                          + I_leak - I_up + I_rel));  
*/			     			    
  // Eq. 96 - TT04
  // from paper: CaSRtotal = Ca_srbufsr - CaSR
  y_TT04[tt04_CaSR] += dt_ms * 
                               (   ((Ca_srbufsr * tt04ps.VC)/tt04ps.VSR) 
			         * (I_up - (I_rel + I_leak))   );

  
  // Eq. 97 - TT04
  // the minus is in the tt04ps.inverseVcF factor
  y_TT04[tt04_Nai]+= (dt_ms * (tt04ps.Cm * tt04ps.inverseVcF * (I_Na + I_bNa + 3 * (I_NaK + I_NaCa))));
  // Eq. 98 - TT04
  // the minus is in the tt04ps.inverseVcF factor
  // in the paper, I_ax (axial current flow) is included but not used in CellML model
  y_TT04[tt04_Ki] += (dt_ms * (tt04ps.Cm * tt04ps.inverseVcF * (I_K1 + I_to + I_Kr + I_Ks - (2 * I_NaK) + I_pK + i_stim)));


     
   
   
   
  // Eq. 3 tenTusscher et al. AJP Heart Circ Physiol 2004
  Iion = I_Na + I_K1 + I_to + I_Kr + I_Ks + I_CaL + I_NaCa + I_NaK + I_pCa + I_pK + I_bCa + I_bNa;
  return(Iion);      
  
  
}

void IBM_TT04::Print()//ostream &tempstr, double tArg,  double V) 
{  
};

void IBM_TT04::LongPrint()//ostream &tempstr, double tArg,  double V) 
{
};
