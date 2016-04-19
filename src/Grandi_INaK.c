#include <math.h>
#include "Grandi.h"
#include "Grandi_INaK.h"

void Grandi_INaKFunc(CELLPARMS *parmsPtr, double *cell, int pOffset, DERIVED *derived, double dt)
{
   VOLTAGE *voltage = (VOLTAGE *)cell; 
   CONCENTRATIONS  *concentrations = (CONCENTRATIONS*) (cell + CONCENTRATIONS_OFFSET); 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double v = voltage->Vm; 
   double Naj = concentrations->Naj; 
   double Nasl= concentrations->Nasl; 
   double Ki  = concentrations->Ki; 

   double KmNaip=11.0*(1.0-0.25*ISO);
   double KmKo =1.5;
   double sigma=(exp(Nao/67.3)-1.0)/7.0;
   double fnak=1.0/(1.0+0.1245*exp(-0.1*v*FRT)+0.0365*sigma*exp(-v*FRT));
   double Fsl=1.0-Fjunc;

   derived->I.NaK_junc=Fjunc*cP->IbarNaK*fnak*Ko /(1.0+pow(KmNaip/Naj,4.0))/(Ko+KmKo);
   derived->I.NaK_sl=Fsl*cP->IbarNaK*fnak*Ko /(1.0+pow(KmNaip/Nasl,4.0)) /(Ko+KmKo);
}
