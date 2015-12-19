// ICab  V,Cai::PCab:
#include <math.h>
#include "OHaraRudy.h"
#include "OHaraRudy_ICab.h"
void OHaraRudy_ICabFunc(CELLPARMS *parmsPtr, double *cell, int pOffset, DERIVED *derived, double dt )
{
   VOLTAGE *voltage = (VOLTAGE *)cell; 
   CONCENTRATIONS  *concentrations = (CONCENTRATIONS*) (cell + CONCENTRATIONS_OFFSET); 
   double V = voltage->Vm; 
   double Cai    =concentrations->Cai;
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double x = zCa*V*FRT; 
   double f;
   xexp(f,x); 
   derived->I.Cab = cP->PCab*zCa*F*(gammaCai*Cai*exp(zCa*V*FRT)-gammaCao*Cao)*f; 
}
