#include "OHaraRudy.h"
#include "OHaraRudy_CaMKtrap.h"
static double aCaMK = 0.05; //1/ms
static double bCaMK = 0.00068; // 1/ms
static double CaMK0 = 0.05; 
static double KmCaMK = 0.15;  // mM
static double KmCaM  = 0.0015;  // mM
void OHaraRudy_CaMKtrapFunc(CELLPARMS *parmsPtr, STATE *state, int pOffset, DERIVED *derived, double dt)
{

   PSTATE *pState = (PSTATE *)(((double *)state)+pOffset) ; 
   double CaMKtrap = pState->CaMKtrap; 
   double Cass = state->Cass; 
   double CaMKBound = CaMK0 * (1.0-CaMKtrap)*Cass/(KmCaM+Cass); 
   double CaMKActive = CaMKBound + CaMKtrap; 
   derived->phiCaMK = CaMKActive/(KmCaMK + CaMKActive);  
   double dCaMKtrap = aCaMK * CaMKBound * CaMKActive - bCaMK*CaMKtrap; 
   pState->CaMKtrap+=dt*dCaMKtrap; 
}
