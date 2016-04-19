#include <math.h>
#include "Grandi.h"
#include "Grandi_IKr.h"

void Grandi_IKrFunc(CELLPARMS *parmsPtr, double *cell, int pOffset, DERIVED *derived, double dt)
{
   VOLTAGE *voltage = (VOLTAGE *)cell; 
   PSTATE *pState = (PSTATE *)(cell+pOffset) ; 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double v = voltage->Vm; 
   double EK = derived->EK; 

   double xkr = pState->xkr; 
   double rkr=1.0/(1.0+exp((v+74.0)/24.0));

   derived->I.Kr=cP->GKr*sqrt(Ko/5.4)*xkr*rkr*(v-EK);

   double xrss=1.0/(1.0+exp(-(v+10.0)/5.0));
   double tauxr=550.0/(1.0+exp((-22.0-v)/9.0))*6.0/(1.0+exp((v-(-11.0))/9.0))+230.0/(1+exp((v-(-40.0))/20.0));

   double dxkr=(xrss-xkr)/tauxr;

   ENDCODE()
     //pState->xkr += dt*dxkr; 
   pState->xkr=xrss-(xrss-xkr)*exp(-dt/tauxr);

}
