// ICab  V,Cai::PCab:
#include <math.h>
#include "OHaraRudy.h"
#include "OHaraRudy_ICab.h"
static double FRT = -1; 
void OHaraRudy_ICabFunc(CELLPARMS *parmsPtr, STATE *state, int pOffset, DERIVED *derived, double dt )
{
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double V = state->Vm; 
   double Cai = state->Cai; 
   derived->I.Cab = cP->PCab*SQ(zCa)*V*FRT*F*(gammaCai*Cai*exp(zCa*V*FRT)-gammaCao*Cao)/(exp(zCa*V*FRT)-1); 
}
COMPONENTINFO OHaraRudy_ICabInit()
{
   if (FRT  < 0) FRT = F/(R*T);
   COMPONENTINFO info;
   info.nVar = nVar; 
   info.varInfo = varInfo;
   info.func = OHaraRudy_ICabFunc;
   info.access = OHaraRudy_ICabAccess;
   return info;
}

