#include <math.h>
#include "Grandi.h"
#include "Grandi_IKs.h"

void Grandi_IKsFunc(CELLPARMS *parmsPtr, double *cell, int pOffset, DERIVED *derived, double dt)
{
   VOLTAGE *voltage = (VOLTAGE *)cell; 
   CONCENTRATIONS  *concentrations = (CONCENTRATIONS*) (cell + CONCENTRATIONS_OFFSET); 

   PSTATE *pState = (PSTATE *)(cell+pOffset) ; 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double v = voltage->Vm; 

   double EKs = derived->EKs; 
   double xks = pState->xks; 

   double phi=(1.0+1.0*AF+2.0*ISO);
   double Fsl=1.0-Fjunc;

   derived->I.Ks_junc=Fjunc*cP->GKs*phi*xks*xks*(v-EKs);
   derived->I.Ks_sl=Fsl*cP->GKs*phi*xks*xks*(v-EKs);

   double xsss=1.0 / (1.0+exp(-(v+40.0*ISO + 3.8)/14.25));
   double tauxs=990.1/(1.0+exp(-(v+40.0*ISO+2.436)/14.12));

   double dxks=(xsss-xks)/tauxs;

   ENDCODE()
     //pState->xks += dt*dxks; 
   pState->xks=xsss-(xsss-xks)*exp(-dt/tauxs);

}
