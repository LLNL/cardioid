// INaL    V:ENa,phiCaMK:GNaL   :mL,hL,hLCaMK
#include <math.h>
#include "OHaraRudy.h"
#include "OHaraRudy_INaL.h"
void OHaraRudy_INaLFunc(CELLPARMS *parmsPtr, STATE *state, int pOffset, DERIVED *derived, double dt)
{
   PSTATE *pState = (PSTATE *)(((double *)state)+pOffset) ; 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double V = state->Vm; 
   double ENa = derived->ENa; 
   double phiCaMK=derived->phiCaMK;
   //   Gates needed to calculate INaL; 
   double mL=pState->mL;
   double hL=pState->hL; 
   double hLCaMK=pState->hLCaMK;     

   derived->I.NaL= cP->GNaL*(V-ENa)*mL*((1-phiCaMK)*hL+phiCaMK*hLCaMK); 

   double mLMhu =  sige(-(V+42.85)/5.264);
   double mTauR = 6.765*exp((V + 11.64)/34.77)+8.552*exp(-(V + 77.42)/5.955); 
   double mLTauR = mTauR; 
   double dmL = (mLMhu-mL)*mLTauR;  // gate

   double hLMhu = sige((V+87.61)/7.488);
   double hLTauR = 0.005; 
   double dhL = (hLMhu-hL)*hLTauR;  // gate

   double hLCaMKMhu = sige((V+93.81)/7.488);
   double hLCaMKTauR = hLTauR/3.0; 
   double dhLCaMK = (hLCaMKMhu-hLCaMK)*hLCaMKTauR;  // gate
   pState->mL += dt*dmL; 
   pState->hL += dt*dhL; 
   pState->hLCaMK += dt*dhLCaMK; 
}
