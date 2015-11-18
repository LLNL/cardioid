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
