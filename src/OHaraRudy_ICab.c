// ICab  V,Cai::PCab:
#include <math.h>
#include "OHaraRudy.h"
#include "OHaraRudy_ICab.h"
void OHaraRudy_ICabFunc(CELLPARMS *parmsPtr, STATE *state, int pOffset, DERIVED *derived, double dt )
{
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double V = state->Vm; 
   double Cai = state->Cai; 
   derived->I.Cab = cP->PCab*SQ(zCa)*V*FRT*F*(gammaCai*Cai*exp(zCa*V*FRT)-gammaCao*Cao)/(exp(zCa*V*FRT)-1); 
}
