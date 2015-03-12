#ifndef OHARARUDY_H
#define OHARARUDY_H
enum CELLTYPES {ENDO_CELL, EPI_CELL, M_CELL};
typedef struct state_st 
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
   double m;           // 1 true
   double hFast;       // 1 true
   double hSlow;       // 1 true
   double j;           // 1 true
   double hCaMKSlow;   // 1 true
   double jCaMK;       // 1 true
   double mL;          // 1 true
   double hL;          // 1 true
   double hLCaMK;      // 1 true
   double a;           // 1 true
   double iFast;       // 1 true
   double iSlow;       // 1 true
   double aCaMK;       // 1 true
   double iCaMKFast;   // 1 true
   double iCaMKSlow;   // 1 true
   double d;           // 1 true
   double fFast;       // 1 true
   double fSlow;       // 1 true
   double fCaFast;     // 1 true
   double fCaSlow;     // 1 true
   double jCa;         // 1 true
   double n;           // 1 true
   double fCaMKFast;   // 1 true
   double fCaCaMKFast; // 1 true
   double XrFast;      // 1 true
   double XrSlow;      // 1 true
   double Xs1;         // 1 true
   double Xs2;         // 1 true
   double XK1;         // 1 true
   double JrelNP;      // mM/ms true
   double JrelCaMK;    // mM/ms true
   double CaMKtrap;    // 1 true
} STATE; 

typedef struct cellParms_st 
{
//        name        units checkpoint
   double GNaFast;       // mS/uF false
   double GNaL;          // mS/uF false
   double Gto;           // mS/uF false
   double PCaL;          // cm/s  false  PCaCaMK  = 1.1*PCaL; 
   double PCaNa;         // cm/s  false  PCaNaCaMK=1.1*PCaNa; 
   double PCaK;          // cm/s  false  PCaKCaMK =1.1*PCaK 
   double GKr ;          // mS/uF false
   double GK1 ;          // mS/uF false
   double GKs ;          // mS/uF false
   double GNaCai;        // uA/uF false
   double GNaCass;       // uA/uF false
   double PNaK;          // mV/mM  false
   double PNab;          // cm/s  false
   double PCab;          // cm/s  false
   double GKb;           // mS/uF false
   double GpCa;          // mS/uF false
   double alphaJrelNP;   // mM/mV false
   double betaJrelNP;    // mM/mS false
   double alphaJrelCaMK; // mM/mV false
   double betaJrelCaMK;  // mM/mS false
   double cJup;          // mM/ms false
   double CMDN;          // mM    false
   double aDelta;        //  1    false
} CELLPARMS;

#ifdef __cplusplus
extern "C" 
{
#endif 
void OHaraRudySetup();
double OHaraRudyIntegrate(double dt, double Istim, STATE *state, CELLPARMS *cP, STATE *D);
double OHaraRudyIntegrateS(double dt, double Istim, STATE *state, CELLPARMS *cP, STATE *D);
void OHaraRudyInitialState(int cellType, STATE *state);
CELLPARMS *OHaraRudyCellParms(int cellType);
#ifdef __cplusplus
}
#endif 
#endif
