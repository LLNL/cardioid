/*
*/
#include "IBM_TT04_LUT.hh"



  char IBM_TT04_LUT::TT04LUT_EVfilename[255];
  TT04ParameterStruct IBM_TT04_LUT::TT04LUT_ttps; //fps stands for TTParameterStruct

  double IBM_TT04_LUT::TT04LUT_Xr1_inf[V_RESOLUTION];
  double IBM_TT04_LUT::TT04LUT_tau_Xr1[V_RESOLUTION];
  double IBM_TT04_LUT::TT04LUT_Xr2_inf[V_RESOLUTION];
  double IBM_TT04_LUT::TT04LUT_tau_Xr2[V_RESOLUTION];
  double IBM_TT04_LUT::TT04LUT_Xs_inf[V_RESOLUTION];
  double IBM_TT04_LUT::TT04LUT_tau_Xs[V_RESOLUTION];

  double IBM_TT04_LUT::TT04LUT_m_inf[V_RESOLUTION];
  double IBM_TT04_LUT::TT04LUT_tau_m[V_RESOLUTION];
  double IBM_TT04_LUT::TT04LUT_h_inf[V_RESOLUTION];
  double IBM_TT04_LUT::TT04LUT_tau_h[V_RESOLUTION];
  double IBM_TT04_LUT::TT04LUT_j_inf[V_RESOLUTION];
  double IBM_TT04_LUT::TT04LUT_tau_j[V_RESOLUTION];
  double IBM_TT04_LUT::TT04LUT_d_inf[V_RESOLUTION];
  double IBM_TT04_LUT::TT04LUT_tau_d[V_RESOLUTION];
  double IBM_TT04_LUT::TT04LUT_f_inf[V_RESOLUTION];
  double IBM_TT04_LUT::TT04LUT_tau_f[V_RESOLUTION];
  double IBM_TT04_LUT::TT04LUT_CaL_P1[V_RESOLUTION];
  double IBM_TT04_LUT::TT04LUT_CaL_P2[V_RESOLUTION];
  double IBM_TT04_LUT::TT04LUT_EpiM_s_inf[V_RESOLUTION];
  double IBM_TT04_LUT::TT04LUT_EpiM_tau_s[V_RESOLUTION];
  double IBM_TT04_LUT::TT04LUT_Endo_s_inf[V_RESOLUTION];
  double IBM_TT04_LUT::TT04LUT_Endo_tau_s[V_RESOLUTION];
  double IBM_TT04_LUT::TT04LUT_r_inf[V_RESOLUTION];
  double IBM_TT04_LUT::TT04LUT_tau_r[V_RESOLUTION];     

    
  // lookup tables
  double IBM_TT04_LUT::TT04LUT_rec_ipK[V_RESOLUTION];
  double IBM_TT04_LUT::TT04LUT_NaCa_P1[V_RESOLUTION];
  double IBM_TT04_LUT::TT04LUT_NaCa_P2[V_RESOLUTION];
  double IBM_TT04_LUT::TT04LUT_NaK_P1[V_RESOLUTION];
                                                                                                                                            
                                                                                                                                                 
     
IBM_TT04_LUT::IBM_TT04_LUT() 
{
};


IBM_TT04_LUT::~IBM_TT04_LUT() 
{ 
};

void IBM_TT04_LUT::TT04LUT_Init()
{

  TT04LUT_InitParameters();
  TT04LUT_InitTable();
           
}


void IBM_TT04_LUT::TT04LUT_InitParameters() 
{

  // constants CellML tenTusscher et al. AJP Heart Circ Physiol 2004
  TT04LUT_ttps.R = 8314.472;         // gas constant                     - J*K^(-1)*mol^(-1)   CellML Note: value different from TT04 publication
  TT04LUT_ttps.T = 310.0;            // temperature                      - K
  TT04LUT_ttps.F = 96485.3415;          // Faraday constant                 - C/mmol           CellML Note: value different from TT04 publication
  TT04LUT_ttps.Cm = 0.185;               // Cell capacitance per unit surface area  - microF/cm^2
  TT04LUT_ttps.VC = 0.016404;           // Cytoplasmic volume               - microm^3         CellML Note: value different from TT04 publication
  TT04LUT_ttps.pKNa = 0.03;          // relative IKs permeability to Na+
  TT04LUT_ttps.Nao = 140.0;            // Extracellular Na+ concentration  - mM
  TT04LUT_ttps.Ko = 5.4;             // Extracellular K+ concentration   - mM
  TT04LUT_ttps.Cao = 2.0;              // Extracellular Ca2+ concentration - mM
  TT04LUT_ttps.GK1 = 5.405;          // Maximal IK1 conductance          - nS/pF
  TT04LUT_ttps.GKr = 0.096;          // Maximal IKr conductance          - nS/pF
  TT04LUT_ttps.GKs_epi_endo = 0.245; // Maximal IKs (epi/endo) conductance - nS/pF
  TT04LUT_ttps.GKs_M = 0.062;        // Maximal IKs (M) conductance        - nS/pF
  TT04LUT_ttps.GKs = 0.0;  // set later to cell type value
  TT04LUT_ttps.GNa = 14.838;         // Maximal INa conductance          - nS/pF
  TT04LUT_ttps.GbNa = 0.00029;       // Maximal IbNa conductance               - nS/pF
  TT04LUT_ttps.GCaL = 0.000175;      // = 1.75*10^(-4) != {0.10662224 = 1.75^(-4)} Maximal ICaL conductance - cm^3*microF^(-1)*s^(-1)
  TT04LUT_ttps.GbCa = 0.000592;      // Maximal IbCa conductance               - nS/pF
  TT04LUT_ttps.Gto_epi_M = 0.294;    // Maximal Gto (epi/M) conductance  - nS/pF
  TT04LUT_ttps.Gto_endo = 0.073;     // Maximal Gto (endo) conductance   - nS/pF
  TT04LUT_ttps.Gto = 0.0;  // will be set later according to cell type
  TT04LUT_ttps.PNaK = 1.362;         // Maximal INaK                           - pA/pF
  TT04LUT_ttps.KmK = 1.0;              // Ko half-saturation constant of INaK    - mM
  TT04LUT_ttps.KmNa = 40.0;            // Nai half-saturation constant of INaK   - mM
  TT04LUT_ttps.kNaCa = 1000;         // Maximal INaCa                          - pA/pF
  TT04LUT_ttps.ksat = 0.1;           // Saturation factor for INaCa
  TT04LUT_ttps.alpha = 2.5;          // Factor enhancing outward nature of INa CA
  TT04LUT_ttps.gamma = 0.35;         // Voltage dependence parameter of INaCa 
  TT04LUT_ttps.KmCa = 1.38;          // Cai half saturation constant for INaCa - mM
  TT04LUT_ttps.KmNai = 87.5;         // Nai half saturation constant for INaCa - mM
  TT04LUT_ttps.GpCa = 0.825;         // Maximal IpCa conductance               - nS/pF    CellML Note: the TT04 publication has wrong number (0.025) cited
  TT04LUT_ttps.KpCa = 0.0005;        // Cai half-saturation constant of IpCa   - mM
  TT04LUT_ttps.GpK = 0.0146;         // Maximal IpK conductance                - nS/pF
  TT04LUT_ttps.arel = 0.016464;        // Maximal CaSR-dependent Irel            - mM/s   CellML Note: value different from TT04 publication
  TT04LUT_ttps.brel = 0.25;          // CaSR half-saturation constant of Irel  - mM
  TT04LUT_ttps.crel = 0.008232;         // Maximal CaSR-independent Irel          - mM/s  CellML Note: value different from TT04 publication
  TT04LUT_ttps.Kup = 0.00025;        // Half-saturation constant of Iup        - mM
  TT04LUT_ttps.Vleak = 0.00008;      // Maximal Ileak                          - ms^(-1)
  TT04LUT_ttps.Vmaxup = 0.000425;    // Maximal Iup                            - mM/ms
  TT04LUT_ttps.Bufc = 0.15;          // Total cytoplasmic buffer concentration - mM
  TT04LUT_ttps.Kbufc = 0.001;        // Cai half-saturation constant for cytoplasmic buffer concentraton - mM
  TT04LUT_ttps.Bufsr = 10.0;           // Total sarcoplasmic buffer concentration - mM
  TT04LUT_ttps.Kbufsr = 0.3;         // CaSR half-saturation constant for sarcoplasmic buffer - mM
  TT04LUT_ttps.VSR = 0.001094;           // Sarcoplasmic reticulum volume    - microm^3  CellML Note: value different from TT04 publication


  TT04LUT_ttps.RToverF = (TT04LUT_ttps.R * TT04LUT_ttps.T )/( TT04LUT_ttps.F );
  TT04LUT_ttps.inverseRToverF = 1 / TT04LUT_ttps.RToverF;
  TT04LUT_ttps.KopKNaNao = TT04LUT_ttps.Ko +( TT04LUT_ttps.pKNa )*( TT04LUT_ttps.Nao );
  TT04LUT_ttps.CaL_4FF_RT = 4 * TT04LUT_ttps.inverseRToverF * TT04LUT_ttps.F;
  TT04LUT_ttps.inverseVcF = (-1.0)*(1./(( TT04LUT_ttps.VC )*( TT04LUT_ttps.F)));
  TT04LUT_ttps.inverse2VcF = 0.5 * TT04LUT_ttps.inverseVcF;
  TT04LUT_ttps.KmNai3 = TT04LUT_ttps.KmNai *( TT04LUT_ttps.KmNai )*( TT04LUT_ttps.KmNai );
  TT04LUT_ttps.Nao3   = TT04LUT_ttps.Nao *( TT04LUT_ttps.Nao )*( TT04LUT_ttps.Nao );

                              
}
                              


void IBM_TT04_LUT::TT04LUT_InitTable()
{

  printf("TableInit\n"); 
  for(double V=-RangeTabhalf+.0001; V<RangeTabhalf; V+=dDivisionTab) 
     {
       //V in mV
       const int Vi=(int)(DivisionTab*(RangeTabhalf+V)+.5);
       
      
       // Fast Na+ current
       // Eq. 28 - TT04
       TT04LUT_m_inf[Vi]=1.0/((1.0+exp((-56.86-V)/9.03))*(1.+exp((-56.86-V)/9.03)));
       // Eq. 29 - TT04
       const double a_m=1.0/(1.0+exp((-60.0-V)/5.0));
       // Eq. 30 - TT04
       const double b_m = (0.1/(1.0+exp((V+35.0)/5.0)))
                        + (0.1/(1.0+exp((V-50.0)/200.0)));
       // Eq. 31 - TT04
       TT04LUT_tau_m[Vi]=a_m*b_m;
       
       // Eq. 32 - TT04
       TT04LUT_h_inf[Vi]=1.0/((1.0+exp((V+71.55)/7.43))*(1.0+exp((V+71.55)/7.43)));
       // Eq. 33 - TT04
       // Eq. 34 - TT04
       double a_h = 0.0;
       double b_h = 0.0;
       if (V >= -40.0)
          {
	    a_h = 0.0;
	    b_h = (0.77/(0.13*(1.+exp((V+10.66)/(-11.1)))));
	  }
       else
          {
	    a_h = (0.057*exp(-(V+80.)/6.8));
	    b_h = (2.7*exp(0.079*V))+(310000*exp(0.3485*V));
	  }
       // Eq. 35 - TT04
       TT04LUT_tau_h[Vi]= 1.0/(a_h + b_h);
	  
      
       
       // Eq. 36 - TT04
       // TT04LUT_j_inf[Vi]=TT04LUT_h_inf[Vi]; // according to old implementation based on IBT code
       TT04LUT_j_inf[Vi] = 1.0/((1+exp((V+71.55)/7.43)) * (1+exp((V+71.55)/7.43)));
       // Eq. 37 - TT04
       // Eq. 38 - TT04
          
       double a_j = 0.0;
       double b_j = 0.0;
       if (V>=-40.)
          {
            a_j=0.0;
            b_j=(0.6*exp(0.057*V)/(1.0+exp(-0.1*(V+32.0))));
          }
       else
          {
            a_j =  (-25428) * exp( 0.2444 *V);     
	    a_j -= (6.948e-6)  * exp(-0.04391*V);
	    a_j *= (V + 37.78);                    
            a_j /= (1.0+exp(0.311*(V+79.23)));     
	    
            b_j=(0.02424*exp(-0.01052*V)/(1.0+exp(-0.1378*(V+40.14))));
          }
       // Eq. 39 - TT04
       TT04LUT_tau_j[Vi]= 1.0/(a_j+b_j);

       
       // L-type Ca2+ current
	// Eq. 40 - TT04
	TT04LUT_CaL_P1[Vi] = TT04LUT_ttps.GCaL * (V * TT04LUT_ttps.CaL_4FF_RT);
	TT04LUT_CaL_P1[Vi] *= exp(2.0 * V * TT04LUT_ttps.inverseRToverF);	
	TT04LUT_CaL_P1[Vi] /= (exp(2.0 * V * TT04LUT_ttps.inverseRToverF)-1.0);
	
	TT04LUT_CaL_P2[Vi] = TT04LUT_ttps.GCaL * (V * TT04LUT_ttps.CaL_4FF_RT);	
        TT04LUT_CaL_P2[Vi] *= (-0.341 * TT04LUT_ttps.Cao);
	TT04LUT_CaL_P1[Vi] /= (exp(2.0 * V * TT04LUT_ttps.inverseRToverF)-1.0);
   
	// Eq. 41 - TT04
        TT04LUT_d_inf[Vi] = 1.0/(1.+exp((-5.0-V)/7.5));
	// Eq. 42 - TT04
        const double a_d = 1.4/(1.0+exp((-35.0-V)/13.0))+0.25;
	// Eq. 43 - TT04
        const double b_d = 1.4/(1.0+exp((V+5.0)/5.0));
	// Eq. 44 - TT04
        const double c_d = 1.0/(1.0+exp((50.0-V)/20.0));
	// Eq. 45 - TT04
        TT04LUT_tau_d[Vi] = (a_d*b_d)+c_d;
	
	// Eq. 46 - TT04
        TT04LUT_f_inf[Vi] = 1.0/(1.0 + exp((V+20.0)/7.0));
	// Eq. 47 - TT04
        TT04LUT_tau_f[Vi]=(1125.0*exp(-((V+27)*(V+27))/240))+80+(165/(1.+exp((25-V)/10)));

       
        // Transient outward current
	// Eq. 56 - TT04
        TT04LUT_r_inf[Vi] = 1./(1.+ exp((20-V)/6.));
	// Eq. 57 - TT04
        TT04LUT_tau_r[Vi] = (9.5 * exp(-((V+40.)*(V+40.))/1800.))+0.8;
	// Eq. 58 - TT04
        TT04LUT_EpiM_s_inf[Vi] = 1.0/(1.0+exp((V+20.0)/5.0));	
	// Eq. 59 - TT04
	TT04LUT_EpiM_tau_s[Vi] = (85.0 * exp(-((V + 45.0) * (V + 45.0))/320.0))
	                       + (5.0/(1.0+exp((V-20.0)/5.0)))
			       + 3.0;
	// Eq. 60 - TT04
        TT04LUT_Endo_s_inf[Vi] = 1.0/(1.0+exp((V+28.0)/5.0));	
	// Eq. 61 - TT04
        TT04LUT_Endo_tau_s[Vi] = (1000.0 * exp(-((V + 67.0) * (V + 67.0))/1000.0)) + 8.0;
		
		
	// Slow delayed rectifier current
	// Eq. 63 - TT04
        TT04LUT_Xs_inf[Vi]=1.0/(1.0+exp((-5.0-V)/14.0));
	// Eq. 64 - TT04
        const double a_Xs=1100.0/(sqrt(1.0+exp((-10.0-V)/6)));
	// Eq. 65 - TT04
        const double b_Xs=1.0/(1.0+exp((V-60.0)/20.0));
	// Eq. 66 - TT04
        TT04LUT_tau_Xs[Vi]=a_Xs*b_Xs;


        // Rapid delayed rectifier current
	// Eq. 68 - TT04
        TT04LUT_Xr1_inf[Vi]=1.0/(1.0+exp((-26.0 - V)/(7.0)));
	// Eq. 69 - TT04
        const double a_Xr1=450.0/(1.0+exp((-45.0 - V)/( 10.0 )));
	// Eq. 70 - TT04
        const double b_Xr1=6.0/(1.0+exp((V + 30.0)/(11.5)));
	// Eq. 71 - TT04
        TT04LUT_tau_Xr1[Vi]=a_Xr1*b_Xr1;
	// Eq. 72 - TT04
        TT04LUT_Xr2_inf[Vi]=1.0/(1.0+exp((V + 88.0)/24.0));
	// Eq. 73 - TT04
        const double a_Xr2=3.0/(1.0+exp((-60.0 - V)/20.0));
	// Eq. 74 - TT04
        const double b_Xr2=1.12/(1.0+exp((V - 60.0)/(20.0)));
	// Eq. 75 - TT04
        TT04LUT_tau_Xr2[Vi]=a_Xr2*b_Xr2;


	
	// Na+/Ca2+ Exchange current
	// Eq. 80 - TT04
        const double NaCaP1 = TT04LUT_ttps.kNaCa 
                            * exp(( TT04LUT_ttps.gamma * V * TT04LUT_ttps.inverseRToverF ))
			    * (TT04LUT_ttps.Cao);
        const double NaCaP2 = TT04LUT_ttps.kNaCa 
                            * exp((( TT04LUT_ttps.gamma - 1) * V * TT04LUT_ttps.inverseRToverF ))
	  		    * ( TT04LUT_ttps.Nao3 ) * TT04LUT_ttps.alpha;
        const double NaCaP3 = 1.0
                            / (  ( TT04LUT_ttps.KmNai3 + TT04LUT_ttps.Nao3 ) 
			       * ( TT04LUT_ttps.KmCa   + TT04LUT_ttps.Cao )
			       * (1.0 + (TT04LUT_ttps.ksat * exp( (TT04LUT_ttps.gamma - 1) * V * TT04LUT_ttps.inverseRToverF )))
			      );
       
        TT04LUT_NaCa_P1[Vi] = NaCaP1 * NaCaP3;
        TT04LUT_NaCa_P2[Vi] = (-1.0) * NaCaP2 * NaCaP3;

        // Na+/K+ pump current
	// Eq. 81 - TT04        
        TT04LUT_NaK_P1[Vi] = (TT04LUT_ttps.PNaK * TT04LUT_ttps.Ko)
	                   / ( (TT04LUT_ttps.Ko +  TT04LUT_ttps.KmK)
			      *(1.0 + (0.1245 * exp(-0.1 * V * TT04LUT_ttps.inverseRToverF)) 
			            + (0.0353 * exp(-1.0 * V * TT04LUT_ttps.inverseRToverF))
				)
			      );	

        TT04LUT_rec_ipK[Vi] = 1.0/(1.0 + exp((25 - V)/5.98));
      
     }// end for
     

}

