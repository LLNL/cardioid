// INaCai   V,Cai,Nai::GNaCai:
// INaCass  V,Cass,Nass::GNaCass:
#include "OHaraRudy.h"
#include "OHaraRudy_INaCaCommon.h"
#include "OHaraRudy_INaCai.h"
void OHaraRudy_INaCaiFunc(CELLPARMS *parmsPtr, STATE *state, int pOffset, DERIVED *derived, double dt)
{
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double V = state->Vm; 
   derived->I.NaCai =  cP->GNaCai *INaCaYFunc(state->Vm,state->Cai, state->Nai   );
}
COMPONENTINFO OHaraRudy_INaCaiInit()
{
   if (FRT < 0) FRT = F/(R*T); 
   COMPONENTINFO info;
   info.nVar = nVar; 
   info.varInfo = varInfo;
   info.func = OHaraRudy_INaCaiFunc;
   info.access = OHaraRudy_INaCaiAccess;
   return info; 
}
