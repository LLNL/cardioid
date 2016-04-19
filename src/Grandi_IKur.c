#include <math.h>
#include "Grandi.h"
#include "Grandi_IKur.h"

void Grandi_IKurFunc(CELLPARMS *parmsPtr, double *cell, int pOffset, DERIVED *derived, double dt)
{
   VOLTAGE *voltage = (VOLTAGE *)cell; 
   CONCENTRATIONS  *concentrations = (CONCENTRATIONS*) (cell + CONCENTRATIONS_OFFSET); 

   PSTATE *pState = (PSTATE *)(cell+pOffset) ; 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double v = voltage->Vm; 

   double EK = derived->EK; 
   double xkur = pState->xkur; 
   double ykur = pState->ykur; 

   double phi=(1.0-0.5*AF)*(1.0+2.0*ISO)*(1.0+0.2*RA);

   derived->I.Kur=cP->GKur*phi*xkur*ykur*(v-EK);

   double xkurss=1.0/(1.0+exp((v+6.0)/-8.6));
   double tauxkur=9.0/(1.0+exp((v+5.0)/12.0))+0.5;
   double ykurss=1.0/(1.0+exp((v+7.5)/10.0));
   double tauykur=590.0/(1.0+exp((v+60.0)/10.0))+3050.0;

   double dxkur=(xkurss-xkur)/tauxkur;
   double dykur=(ykurss-ykur)/tauykur;

   ENDCODE()
     //pState->xkur += dt*dxkur; 
     //pState->ykur += dt*dykur; 
   pState->xkur=xkurss-(xkurss-xkur)*exp(-dt/tauxkur);
   pState->ykur=ykurss-(ykurss-ykur)*exp(-dt/tauykur);

}
