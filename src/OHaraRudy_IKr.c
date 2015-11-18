// IKr    V:EK:GKr:XrSlow,XrFast
#include <math.h>
#include "OHaraRudy.h"
#include "OHaraRudy_IKr.h"

void OHaraRudy_IKrFunc(CELLPARMS *parmsPtr, STATE *state, int pOffset, DERIVED *derived, double dt)
{

   PSTATE *pState = (PSTATE *)(((double *)state)+pOffset) ; 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double V = state->Vm; 
   double EK = derived->EK; 
   double XrFast = pState->XrFast; 
   double XrSlow = pState->XrSlow; 
   double AXrFast = sige((V+54.81)/38.21); 
   double AXrSlow = 1.0 - AXrFast; 
   double Xr = AXrFast*XrFast + AXrSlow*XrSlow; 
   double RKr = sige((V+55)/75)*sige((V-10)/30); 
   derived->I.Kr = cP->GKr * (V-EK) *sqrt(Ko/5.4) * Xr * RKr; 

   double XrMhu = sige(-(V+8.337)/6.789); 
   double XrFastTau = 12.98  + 1/(0.3652*exp((V-31.66)/3.869) + 4.123e-5*exp(-(V-47.78)/20.38));
   double XrSlowTau =  1.865 + 1/(0.06629*exp((V-34.70)/7.355) + 1.128e-5*exp(-(V-29.74)/25.94));
   double XrFastTauR = 1/XrFastTau; 
   double XrSlowTauR = 1/XrSlowTau; 
   double dXrSlow = (XrMhu - XrSlow)*XrSlowTauR;  // gate
   double dXrFast = (XrMhu - XrFast)*XrFastTauR;  // gate
   pState->XrSlow += dt*dXrSlow; 
   pState->XrFast += dt*dXrFast; 
}
