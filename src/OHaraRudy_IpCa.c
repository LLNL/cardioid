// IpCa = Cai::GpCa:
#include <math.h>
#include "OHaraRudy.h"
#include "OHaraRudy_IpCa.h"

void OHaraRudy_IpCaFunc(CELLPARMS *parmsPtr, STATE *state, int pOffset, DERIVED *derived, double dt )
{
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double Cai = state->Cai; 
   derived->I.pCa = cP->GpCa * Cai/(0.0005 + Cai); 
}
COMPONENTINFO OHaraRudy_IpCaInit()
{
   COMPONENTINFO info;
   info.nVar = nVar; 
   info.varInfo = varInfo;
   info.func = OHaraRudy_IpCaFunc;
   info.access = OHaraRudy_IpCaAccess;
   return info;
}

