#include "OHaraRudy.h"
#include "OHaraRudy_Fluxes.h"
static double tauRDiffNa=0.5; // ms
static double tauRDiffK =0.5; // ms
static double tauRDiffCa=5.0;// ms
static double tauRtr    = 0.01;// ms
static double deltaKmPLB = 0.00017 ; // mM
static double deltaJupCaMK=1.75; 
void OHaraRudy_FluxesFunc(CELLPARMS *parmsPtr, STATE *state, int pOffset, DERIVED *derived, double dt )
{
   FLUXES J; 
   PSTATE *pState = (PSTATE *)(((double *)state)+pOffset) ; 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double phiRelCaMK=derived->phiCaMK;
   double Cajsr = state->Cajsr; 
   double Cansr = state->Cansr; 
   double Nass = state->Nass; 
   double Cass = state->Cass; 
   double Cai = state->Cai; 
   double Nai = state->Nai; 
   double Ki = state->Ki; 
   double Kss = state->Kss; 
   derived->J.rel = (1.0-phiRelCaMK)*pState->JrelNP + phiRelCaMK*pState->JrelCaMK; 

   double x = SQ(SQ(SQ(Cajsr/1.5))); 
   double y = (Cajsr/0.123); 
   double alpha0   = (-derived->I.CaL)*x/(1.0+x);
   double beta0    = y/(1+y); 

   double mhuJrelNP  = cP->alphaJrelNP  *alpha0 ;
   double tauJrelNP  = cP->betaJrelNP   *beta0;
   tauJrelNP = MAX(tauJrelNP,0.001); 
   double dJrelNP   = (mhuJrelNP  -pState->JrelNP  )/tauJrelNP; 

   double mhuJrelCaMK= cP->alphaJrelCaMK*alpha0;
   double tauJrelCaMK= cP->betaJrelCaMK *beta0; 
   tauJrelCaMK = MAX(tauJrelCaMK,0.001); 
   double dJrelCaMK = (mhuJrelCaMK-pState->JrelCaMK)/tauJrelCaMK; 

   double JupNP = cP->cJup*Cai/(0.00092 + Cai) ; 

   double JupCaMK = cP->cJup*(1 + deltaJupCaMK)*Cai/(0.00092 -deltaKmPLB + Cai) ;

   double phiUpCaMK=phiRelCaMK;
   double Jleak = 0.0039375*Cansr/15.0;
   derived->J.up = (1-phiUpCaMK)*JupNP + phiUpCaMK*JupCaMK - Jleak; 

   derived->J.diffNa = (Nass-Nai   )*tauRDiffNa;
   derived->J.diffCa = (Cass-Cai   )*tauRDiffCa;
   derived->J.diffK  = (Kss-Ki     )*tauRDiffK;
   derived->J.tr     = (Cansr-Cajsr)*tauRtr; 

   pState->JrelNP  += dt*dJrelNP; 
   pState->JrelCaMK+= dt*dJrelCaMK; 

}
