// INaCai   V,Cai,Nai::GNaCai:
// INaCass  V,Cass,Nass::GNaCass:
#include "OHaraRudy.h"
#include "OHaraRudy_INaCaCommon.h"
#include "OHaraRudy_INaCass.h"


void OHaraRudy_INaCassFunc(CELLPARMS *parmsPtr, STATE *state, int pOffset, DERIVED *derived, double dt)
{
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double V = state->Vm; 
   derived->I.NaCass = cP->GNaCass*INaCaYFunc(state->Vm,state->Cass, state->Nass  );
}

COMPONENTINFO OHaraRudy_INaCassInit()
{
   if (FRT < 0) FRT = F/(R*T); 
   COMPONENTINFO info;
   info.nVar = nVar; 
   info.varInfo = varInfo;
   info.func = OHaraRudy_INaCassFunc;
   info.access = OHaraRudy_INaCassAccess;
   return info;
}
