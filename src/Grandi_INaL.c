#include <math.h>
#include "Grandi.h"
#include "Grandi_INaL.h"
void Grandi_INaLFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt)
{
   VOLTAGE *voltage = (VOLTAGE *)state; 
   PSTATE *pState = (PSTATE *)(state+pOffset) ; 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double v = voltage->Vm; 
   double ENa_junc = derived->ENa_junc; 
   double ENa_sl = derived->ENa_sl; 
   double Fsl=1.0-Fjunc;

   double mL=pState->mL;
   double hL=pState->hL; 

   derived->I.NaL_junc=Fjunc*cP->GNaL*mL*mL*mL*hL*(v-ENa_junc);
   derived->I.NaL_sl=Fsl*cP->GNaL*mL*mL*mL*hL*(v-ENa_sl);

   double aml=0.32*(v+47.13)/(1.0-exp(-0.1*(v+47.13)));
   double bml=0.08*exp(-v/11.0);
   double hlinf=1.0/(1.0+exp((v+91.0)/6.1));
   double tauhl=600.0;

   double dmL=aml*(1.0-mL)-bml*mL;
   double dhL=(hlinf-hL)/tauhl;

   ENDCODE()
     //pState->mL += dt*dmL; 
     //pState->hL += dt*dhL; 

   double taumL=1.0/(aml+bml);
   double mLss=aml*taumL;
   pState->mL=mLss-(mLss-mL)*exp(-dt/taumL);
   pState->hL=hlinf-(hlinf-hL)*exp(-dt/tauhl);

}
