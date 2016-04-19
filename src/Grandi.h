#ifndef GRANDI_H
#define GRANDI_H
#include <stdio.h>
#include <math.h>
#define doEnd

#define MAX(A, B) ((A) > (B) ? (A) : (B))
#define sigm(x)   (1.0/(1.0 + (x)))
#define sigm2(x)  (1.0/SQ(1.0 + (x)))
#define sige(x)   (1.0/(1.0 + exp((x))))
#define SQ(x)   ((x)*(x))
#define CUBE(x)   ((x)*(x)*(x))

static double AF=0.0;
static double ISO=0.0;
static double RA=0.0;

//Physical Constants 
//
static double R=8314.0;       // [J/kmol*K]
static double F=96485.0;   // [C/mol]
static double T=310.0;     // [K]
static double FRT=-1;

// Cell geometry
static double cellLength=100.0;     // cell length [um]113;//100
static double cellRadius=10.25;   // cell radius [um]12;//10.25

// Fractional currents in compartments
static double Fjunc=0.11;
static double Fjunc_CaL=0.9;

static double pNaK = 0.01833; 
static double Cli=15.0;   // Intracellular Cl  [mM]
static double Clo=150.0;  // Extracellular Cl  [mM]
static double Ko=5.4;   // Extracellular K   [mM]
static double Nao=140.0;  // Extracellular Na  [mM]
static double Cao=1.8;  // Extracellular Ca  [mM]
static double Mgi=1.0;    // Intracellular Mg  [mM]

static double Cmem=1.1e-10;   // [F] membrane capacitance 1.3810e-10;//

static double Q10SRCaP=2.6;          // [none]
static double Vmax_SRCaP=5.3114e-3;  // [mM/msec] (286 umol/L cytosol/sec)
static double Kmr=1.7;               // [mM]L cytosol
static double hillSRCaP=1.787;       // [mM]
static double ks=25.0;                 // [1/ms]
static double kom=0.06;              // [1/ms]
static double kiCa=0.5;              // [1/mM/ms]
static double kim=0.005;             // [1/ms]
static double ec50SR=0.45;           // [mM]
static double MaxSR=15.0;
static double MinSR=1.0;
static double Bmax_TnClow=70.0e-3;    // [mM]                      // TnC low affinity
static double kon_tncl=32.7;        // [1/mM/ms]
static double Bmax_TnChigh=140e-3;  // [mM]                      // TnC high affinity
static double koff_tnchca=0.032e-3; // [1/ms]
static double kon_tnchca=2.37;      // [1/mM/ms]
static double koff_tnchmg=3.33e-3;  // [1/ms]
static double kon_tnchmg=3.0e-3;      // [1/mM/ms]
static double Bmax_CaM=24.0e-3;       // [mM] **? about setting to 0 in c-code**   // CaM buffering
static double koff_cam=238.0e-3;      // [1/ms]
static double kon_cam=34.0;           // [1/mM/ms]
static double Bmax_myosin=140.0e-3;   // [mM]                      // Myosin buffering
static double koff_myoca=0.46e-3;   // [1/ms]
static double kon_myoca=13.8;       // [1/mM/ms]
static double koff_myomg=0.057e-3;  // [1/ms]
static double kon_myomg=0.0157;     // [1/mM/ms]
static double koff_sr=60.0e-3;        // [1/ms]
static double kon_sr=100.0;           // [1/mM/ms]
static double koff_sll=1300.0e-3;     // [1/ms]
static double kon_sll=100.0;          // [1/mM/ms]
static double koff_slh=30.0e-3;       // [1/ms]
static double kon_slh=100.0;          // [1/mM/ms]
static double koff_csqn=65.0;         // [1/ms]
static double kon_csqn=100.0;         // [1/mM/ms]
static double Bmax_Naj=7.561;       // [mM] // Na buffering
static double Bmax_Nasl=1.65;       // [mM]
static double koff_na=1e-3;         // [1/ms]
static double kon_na=0.1e-3;        // [1/mM/ms]

#define CELLPARMS double

enum CELLTYPES {RA_SR, LA_SR, RA_AF, LA_AF};
enum varTypes { PARAMETER_TYPE,PSTATE_TYPE,END_VARINFO}; 
enum accessTypes { READ, WRITE}; 
typedef struct voltage_st 
{
   double Vm;          // mV
   double dVm;         // mV/msec
   double iStim;       // mV/msec
} VOLTAGE; 
#define CONCENTRATIONS_OFFSET 3
typedef struct concentrations_st 
{ 
//        name         units checkpoint
   double Csqnb;
   double NaBj;
   double NaBsl;
   double Naj;         // mM true
   double Nasl;         // mM true
   double Nai;         // mM true
   double Ki;         // mM true
   double Casr;         // mM true
   double Caj;         // mM true
   double Casl;         // mM true
   double Cai;         // mM true
}  CONCENTRATIONS; 

typedef struct flux_st    { double rel, up, leak, CaB_cytosol, CaB_junc, CaB_sl;} FLUXES; 
typedef struct current_st 
{ 
  double stimulus;   
  double Na_junc;
  double Na_sl;
  double NaL_junc;
  double NaL_sl;
  double Nab_junc;
  double Nab_sl;
  double NaK_junc;
  double NaK_sl;
  double Kr;
  double Ks_junc;
  double Ks_sl;
  double Kur;
  double Kp_junc;
  double Kp_sl;
  double to;
  double K1;
  double ClCa_junc;
  double ClCa_sl;
  double Clb;
  double Ca_junc;
  double Ca_sl;
  double Cab_junc;
  double Cab_sl;
  double CaNa_junc;
  double CaNa_sl;
  double CaK;
  double NCX_junc;
  double NCX_sl;
  double pCa_junc;
  double pCa_sl;

} CURRENTS;

typedef struct derived_st 
{ 
  double ENa_junc, ENa_sl, EK, EKs, ECa_junc, ECa_sl, ECl; 
  CURRENTS I; 
  FLUXES J;
  double *dState;
  double Vo; 
} DERIVED; 

typedef struct varInfo_str 
{
   char *name; 
   int type;
   int index;
   double defaultValueRA_SR;
   double defaultValueLA_SR;
   double defaultValueRA_AF;
   double defaultValueLA_AF;
   char *units;
} VARINFO; 

typedef struct componentInfo_st
{
   char *compName; 
   int pStateSize; 
   int parmsSize; 
   int nVar; 
   VARINFO *varInfo;
   void (*func)(CELLPARMS *, double *, int, DERIVED *, double); 
   void (*access)(int, int, double *, double *, double *); 
} COMPONENTINFO;
#define xexp(f,x) do {if ((x)*(x) > 1e-10)  f = (x)/(exp((x))-1.0);  else   f = 1/(1+0.5*(x)+ 0.1666666666666667*(x)*(x));  } while(0)
#ifdef __cplusplus
extern "C" 
{
#endif 
void GrandiInit(double dt, int nCells, int *cellType);
int GrandiGet_nComp(); 
COMPONENTINFO* GrandiGet_compInfo(); 
void GrandiSetValue(int, int handle, double cell); 
double GrandiGetValue(int, int handle); 
void GrandiCalc(); 
void  GrandiGet(double *dVm);
void  GrandiPut(double dt, int nCells, const double *Vm, double const *iStim);
#ifdef __cplusplus
}
#endif 
#endif
