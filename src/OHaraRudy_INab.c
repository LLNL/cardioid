// INab  V,Nai:K:PNab:
#include <math.h>
#include "OHaraRudy.h"
#include "OHaraRudy_INab.h"

static double FRT = -1; 

void OHaraRudy_INabFunc(CELLPARMS *parmsPtr, STATE *state, int pOffset, DERIVED *derived,  double dt)
{
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double V = state->Vm; 
   double Nai = state->Nai; 
   derived->I.Nab = cP->PNab*SQ(zNa) *V*FRT*F*(Nai*exp(zNa*V*FRT)-Nao)/(exp(zNa*V*FRT)-1); 
}
COMPONENTINFO OHaraRudy_INabInit()
{
   if (FRT  < 0) FRT = F/(R*T);
   COMPONENTINFO info;
   info.nVar = nVar; 
   info.varInfo = varInfo;
   info.func = OHaraRudy_INabFunc;
   info.access = OHaraRudy_INabAccess;
   return info;
}

