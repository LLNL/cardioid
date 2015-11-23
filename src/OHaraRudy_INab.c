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
   derived->I.Nab = cP->PNab*SQ(zNa) *V*FRT*F*(Nai*exp(zNa*V*FRT)-Nao)/(exp(zNa*V*FRT)-1); 
}
