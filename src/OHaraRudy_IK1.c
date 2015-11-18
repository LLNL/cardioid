// IK1  V:EK:GK1:XK1
#include <math.h>
#include "OHaraRudy.h"
#include "OHaraRudy_IK1.h"
void OHaraRudy_IK1Func(CELLPARMS *parmsPtr, STATE *state, int pOffset, DERIVED *derived, double dt)
{

   PSTATE *pState = (PSTATE *)(((double *)state)+pOffset) ; 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double V = state->Vm; 
   double EK = derived->EK; 
   double XK1 = pState->XK1; 

   double RK1 = sige((V+105.8-2.6*Ko)/9.493); 
   double phi = sqrt(Ko)*RK1; 
   derived->I.K1 = cP->GK1 * (V-EK) * XK1 * phi;

   double XK1Mhu = sige(-(V + 2.5538*Ko + 144.59)/(1.5692*Ko + 3.8115)); 
   double XK1TauR = (exp(-(V+127.2)/20.36) + exp((V+236.8)/69.33))/122.2;
   double dXK1 = (XK1Mhu-XK1)*XK1TauR; // gate
   pState->XK1 += dt*dXK1; 
}
