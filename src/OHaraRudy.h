#ifndef OHARARUDY_H
#define OHARARUDY_H
#include <stdio.h>

#define MAX(A, B) ((A) > (B) ? (A) : (B))
#define sigm(x)   (1.0/(1.0 + (x)))
#define sigm2(x)  (1.0/SQ(1.0 + (x)))
#define sige(x)   (1.0/(1.0 + exp((x))))
#define SQ(x)   ((x)*(x))
#define CUBE(x)   ((x)*(x)*(x))

//Physical Constants 
//
static double R   = 8314.0;  //J/kmol/K
static double F   = 96485; //Coulomb/mol

static double T   = 310; //K
static double FRT = -1; 
static double Cm     = 1.0    ; // uF; 


static double PRNaK = 0.01833; 
static double Nao =  140; //mM;
static double Cao =  1.8; // mM
static double Ko  =  5.4; // mM

static double gammaCai= 1.0 ;
static double gammaCao= 0.341; 
static double zCa = 2;
static double gammaNai= 0.75;
static double gammaNao= 0.75; 
static double zNa = 1;
static double gammaKi= 0.75;
static double gammaKo= 0.75; 
static double zK = 1;
static double AfFast=0.6; 

#define CELLPARMS double

enum CELLTYPES {ENDO_CELL, EPI_CELL, M_CELL};
enum varTypes { PARAMETER_TYPE,PSTATE_TYPE,END_VARINFO}; 
enum accessTypes { READ, WRITE}; 
typedef struct commonState_st 
{ 
//        name         units checkpoint
   double Vm;          // mV false
   double Nai;         // mM true
   double Nass;        // mM true
   double Ki;          // mM true
   double Kss;         // mM true
   double Cai;         // mM true
   double Cass;        // mM true
   double Cansr;       // mM true
   double Cajsr;       // mM true
}  STATE; 

typedef struct flux_st    { double diffNa, diffK, diffCa, rel, up, tr;} FLUXES; 
typedef struct current_st 
{ 
   double stimulus; 
   double NaCai;
   double NaCass;
   double NaK;
   double Nab;
   double Cab;
   double Kb;
   double pCa;

   double NaFast; 
   double NaL; 
   double to;
   double Kr;
   double Ks;
   double K1;
   double CaL;
   double CaNa;
   double CaK;

} CURRENTS;

typedef struct derived_st 
{ 
  double dVm; 
  double ENa, EK, EKs; 
  double phiCaMK; 
  CURRENTS I; 
  FLUXES J;
} DERIVED; 

typedef struct varInfo_str 
{
   char *name; 
   int type;
   int index;
   double defaultValueENDO;
   double defaultValueEPI;
   double defaultValueM;
   char *units;
} VARINFO; 

typedef struct componentInfo_st
{
   int pStateSize; 
   int parmsSize; 
   int nVar; 
   VARINFO *varInfo;
   void (*func)(CELLPARMS *, STATE *, int, DERIVED *, double); 
   void (*access)(int, int, double *, double *, double *); 
} COMPONENTINFO;
#ifdef __cplusplus
extern "C" 
{
#endif 
void OHaraRudyInit(); 
int OHaraRudyGet_nComp(); 
COMPONENTINFO* OHaraRudyGet_compInfo(); 
void OHaraRudySetValue(int, int handle, double cell); 
void OHaraRudyCalc(); 
void  OHaraRudyGet(double *dVm);
void  OHaraRudyPut(int nCells, const double *Vm, double const *iStim);
#ifdef __cplusplus
}
#endif 
#endif
