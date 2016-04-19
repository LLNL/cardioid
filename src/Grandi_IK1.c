#include <math.h>
#include "Grandi.h"
#include "Grandi_IK1.h"
void Grandi_IK1Func(CELLPARMS *parmsPtr, double *cell, int pOffset, DERIVED *derived, double dt)
{
   VOLTAGE *voltage = (VOLTAGE *)cell; 

   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 

   double v = voltage->Vm; 
   double EK = derived->EK; 

   double phi=(1.0+1.0*AF);
   double aki=1.02/(1+exp(0.2385*(v-EK-59.215)));
   double bki =(0.49124*exp(0.08032*(v+5.476-EK))+exp(0.06175*(v-EK-594.31)))/(1.0+exp(-0.5143*(v-EK+4.753)));

   double kiss=aki/(aki+bki);

   derived->I.K1=cP->GK1*phi*sqrt(Ko/5.4)*kiss*(v-EK);
}
