// INab  V,Nai:K:PNab:
#include <math.h>
#include "OHaraRudy.h"
#include "OHaraRudy_INab.h"

void OHaraRudy_INabFunc(CELLPARMS *parmsPtr, double *cell, int pOffset, DERIVED *derived,  double dt)
{
   VOLTAGE *voltage = (VOLTAGE *)cell; 
   CONCENTRATIONS   *concentrations = (CONCENTRATIONS*) (cell + CONCENTRATIONS_OFFSET); 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double V = voltage->Vm; 
   double Nai = concentrations->Nai; 
   double x = zNa*V*FRT; 
   double f;
   xexp(f,x) ;
   derived->I.Nab = cP->PNab*zNa*F*(Nai*exp(x)-Nao)*f; 
}
