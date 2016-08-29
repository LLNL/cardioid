#include <math.h>
#include "Grandi.h"
#include "Grandi_Ito.h"

void Grandi_ItoFunc(CELLPARMS *parmsPtr, double *cell, int pOffset, DERIVED *derived, double dt )
{
   VOLTAGE *voltage = (VOLTAGE *)cell; 
   PSTATE *pState = (PSTATE *)(((double *)cell)+pOffset) ; 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double v = voltage->Vm; 
   double EK = derived->EK; 

   double xtof = pState->xtof;
   double ytof = pState->ytof; 
   double phi=(1.0-0.7*cP->AF);

   derived->I.to=phi*cP->Gto*xtof*ytof*(v-EK);

   double xtofss=1.0/(1.0+exp(-(v+1.0)/11.0 ));
   double tauxtof=3.5*exp(-(pow(v/30.0,2.0)))+1.5;
   double ytofss=1.0/(1.0+exp((v+40.5)/11.5));
   double tauytof=25.635*exp(-(pow((v+52.45)/15.8827,2.0)))+24.14;

   double dxtof=(xtofss-xtof)/tauxtof;
   double dytof=(ytofss-ytof)/tauytof;

   ENDCODE()
     //pState->xtof += dt*dxtof; 
     //pState->ytof += dt*dytof;
   pState->xtof=xtofss-(xtofss-xtof)*exp(-dt/tauxtof);
   pState->ytof=ytofss-(ytofss-ytof)*exp(-dt/tauytof);

}
