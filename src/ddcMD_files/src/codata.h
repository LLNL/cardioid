/* $Id$  */

#include <math.h>

#ifdef __cplusplus
extern "C"{
#endif

#ifndef CODATA_H
#define CODATA_H
#define CODATA_2006 
//#define OLDDATA 
/* NIST CODATA 2002 */
#ifdef CODATA_2002 
#define me_MKS 9.1093826e-31  /* Electron Mass Standard Uncertainty 0.000 0016e-31 kg CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?me */
#define u_MKS 1.66053886e-27  /*  unified atomic mass unit, Standard Uncertainty 0.000 000 28e-27 kg CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?u */
#define a0_MKS 0.5291772108e-10 /* Bohr radius Standard Uncertainty 0.000 000 0018e-10 meters CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0 */
#define Rinfhc_MKS 2.17987209e-18 /* Standard Uncertainty 0.000 000 37e-18 j CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0 */
#define Ry_eV 13.6056923          /* Standard Uncertainty  0.0000012 eV CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?rydhcev */
#define Eh_eV 27.2113845    /* Standard Uncertainty  0.0000023 eV  CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?threv */
#define hbar_over_Eh_MKS  2.418884326505e-17  /* Standard Uncertainty  0.000000000016e-17 Sec   CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?aut */
#define Na  6.0221415e23
#endif 

/* NIST CODATA 2006 */
#ifdef CODATA_2006 
#define re_MKS     2.8179402894e-15  /* classical electron radius, Standard Uncertainty 0.000 000 0058e-15 m http://physics.nist.gov/cgi-bin/cuu/Value?re */
#define mp_MKS     1.672621777e-27  /*  proton mass, Standard Uncertainty 0.000 000 074 x 10-27 kg http://physics.nist.gov/cgi-bin/cuu/Value?mp */
#define mn_MKS     1.674927351e-27  /*  proton mass, Standard Uncertainty 0.000 000 074 x 10-27 kg http://physics.nist.gov/cgi-bin/cuu/Value?mn */
#define mD_MKS     3.34358348e-27  /*  du mass, Standard Uncertainty 0.000 000 15 x 10-27 kg http://physics.nist.gov/cgi-bin/cuu/Value?md */
#define mT_MKS     5.00735630e-27  /*  triton mass, Standard Uncertainty 0.000 000 22 x 10-27 kg http://physics.nist.gov/cgi-bin/cuu/Value?mt */
#define mh_MKS     5.00641234e-27  /*  He3 mass, Standard Uncertainty 0.000 000 22 x 10-27 kg http://physics.nist.gov/cgi-bin/cuu/Value?mh */
#define mAlpha_MKS 6.64465675e-27  /*  alpha mass, Standard Uncertainty 0.000 000 29 x 10-27 kg http://physics.nist.gov/cgi-bin/cuu/Value?mal */
#define me_MKS     9.10938215e-31  /* Electron Mass, Standard Uncertainty 0.000 000 45e-31 kg CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?me */
#define mec2_MeV   0.510998910   /* electron mass in MeV  0.000 000 013 MeV http://physics.nist.gov/cgi-bin/cuu/CCValue?mec2mev */
#define mec2_MKS   8.18710438e-14   /* electron mass   8.187 104 38e-14h ttp://physics.nist.gov/cgi-bin/cuu/CCValue?mec2 */
#define u_MKS 1.660538782e-27  /*  unified atomic mass unit, Standard Uncertainty 0.000 000 083e-27 kg CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?u */
#define a0_MKS 0.52917720859e-10 /* Bohr radius, Standard Uncertainty 0.000 000 000 36e-10 meters CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0 */
#define Rinfhc_MKS 2.17987197e-18 /*  Rydberg constant times hc in Joules,  Standard Uncertainty 0.000 000 11e-18 j CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0 */
#define Rinfhc_eV 13.60569193       /*  Rydberg constant times hc in eV, Standard Uncertainty  0.000 000 34 eV CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?rydhcev */
#define Eh_eV 27.21138386    /* Hartree Standard Uncertainty  0.000 000 68 eV  CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?threv */
#define Eh 4.35974394e-18    /* Hartree Standard Uncertainty  0.000 000 22e-18 J   CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?hr */
#define hbar_over_Eh_MKS   2.418884326505e-17  /*  atomic unit of time hbar/Eh ,Standard Uncertainty  0.000000000016e-17 Sec   CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?aut */
#define hbar_MKS 1.054571628e-34 /* Planck constant over 2 pi Standard Uncertainty 0.000 000 053e-34 J s CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?hbar */ 
#define NA  6.02214179e23    /*  Avogadro constant Standard Uncertainty  0.00000030e23 Sec   CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?na */
#define Rinf  10973731.568527   /* Rydberg constant  Standard Uncertainty  0.000 073 m^-1  CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?ryd */
#define kB_eV 8.617343e-5  /*  Boltzmann constant in eV/K  Standard Uncertainty 0.000 015 x 10-5 eV/K CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?tkev */
#define kB_MKS 1.3806504e-23  /*  Boltzmann constant in MKS  Standard Uncertainty 0.000 0024 x 10-23 J/K CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?k */
#define AlphaFineStructure  7.2973525376e-3 /* Fine-structure constant Standard Uncertainty 0.000 000 0050e-3 http://physics.nist.gov/cgi-bin/cuu/Value?alph */ 
#define e_MKS 1.602176487e-19  /*  elementary charge Standard Uncertainty 0.000 000 040 x 10-19 C  CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?e */
#define c_MKS 299792458.0  /*  speed of light in vacumm  (exact) m/s  CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?c */
#define mu0_MKS  ((4e-7)*M_PI)  /* magnetic constant (exact) N/A^-2 CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?mu0 */
#define ke_MKS   (c_MKS*c_MKS*(1e-7)) /* Coulomb force constant  mu0*c^2/ (4*Pi) */
#define h_MKS  6.62606896e-34  /* Planck Constant Standard uncertainty 	 0.000 000 33e-34 J s http://physics.nist.gov/cgi-bin/cuu/Value?h */ 
#endif 

#ifdef OLDDATA 
#define re_MKS 2.8179402894e-15  /* classical electron radius, Standard Uncertainty 0.000 000 0058e-15 m http://physics.nist.gov/cgi-bin/cuu/Value?re */
#define me_MKS 9.10938215e-31  /* Electron Mass, Standard Uncertainty 0.000 000 45e-31 kg CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?me */
#define mec2_MeV 0.510998910   /* electron mass in MeV  0.000 000 013 MeV http://physics.nist.gov/cgi-bin/cuu/CCValue?mec2mev */
#define mec2_MKS 8.18710438e-14   /* electron mass   8.187 104 38e-14h ttp://physics.nist.gov/cgi-bin/cuu/CCValue?mec2 */
#define u_MKS 1.660538782e-27  /*  unified atomic mass unit, Standard Uncertainty 0.000 000 083e-27 kg CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?u */
#define a0_MKS 0.52917e-10 /* OLD */
#define Rinfhc_MKS 2.17987197e-18 /*  Rydberg constant times hc in Joules,  Standard Uncertainty 0.000 000 11e-18 j CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0 */
#define Rinfhc_eV 13.60568       /* OLD */
#define Eh_eV 27.21138386    /* Hartree Standard Uncertainty  0.000 000 68 eV  CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?threv */
#define Eh 4.35974394e-18    /* Hartree Standard Uncertainty  0.000 000 22e-18 J   CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?hr */
#define hbar_over_Eh_MKS   2.418884326505e-17  /*  atomic unit of time hbar/Eh ,Standard Uncertainty  0.000000000016e-17 Sec   CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?aut */
#define hbar_MKS 1.054571628e-34 /* Planck constant over 2 pi Standard Uncertainty 0.000 000 053e-34 J s CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?hbar */ 
#define NA  6.02214179e23    /*  Avogadro constant Standard Uncertainty  0.00000030e23 Sec   CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?na */
#define Rinf  10973731.568527   /* Rydberg constant  Standard Uncertainty  0.000 073 m^-1  CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?ryd */
#define kB_eV 8.617343e-5  /*  Boltzmann constant in eV/K  Standard Uncertainty 0.000 015 x 10-5 eV/K CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?tkev */
#define kB_MKS 1.3806504e-23  /*  Boltzmann constant in MKS  Standard Uncertainty 0.000 0024 x 10-23 J/K CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?k */
#define AlphaFineStructure  7.2973525376e-3 /* Fine-structure constant Standard Uncertainty 0.000 000 0050e-3 http://physics.nist.gov/cgi-bin/cuu/Value?alph */ 
#define e_MKS 1.602176487e-19  /*  elementary charge Standard Uncertainty 0.000 000 040 x 10-19 C  CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?e */
#define c_MKS 299792458.0  /*  speed of light in vacumm  (exact) m/s  CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?c */
#define mu0_MKS  ((4e-7)*M_PI)  /* magnetic constant (exact) N/A^-2 CODATA NIST http://physics.nist.gov/cgi-bin/cuu/Value?mu0 */
#define ke_MKS   (c_MKS*c_MKS*(1e-7)) /* Coulomb force constant  mu0*c^2/ (4*Pi) */
#define h_MKS  6.62606896e-34  /* Planck Constant Standard uncertainty 	 0.000 000 33e-34 J s http://physics.nist.gov/cgi-bin/cuu/Value?h */ 
#endif 

extern double kB; 
extern double hbar; 
extern double mu0; 
extern double ke;
extern double hc; 
extern double mec2; 
extern double re; 
extern double qelectron; 
extern double M_e; 
extern double M_p; 
extern double M_n; 
extern double M_D; 
extern double M_T; 
extern double M_He3; 
extern double M_He4 ;
#ifdef __cplusplus
}
#endif

#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
