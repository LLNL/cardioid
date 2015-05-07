// IKs  V,Cai:EKs:GKs:Xs1,Xs2
#include <math.h>
#include "OHaraRudy.h"
#include "OHaraRudy_IKs.h"

void OHaraRudy_IKsFunc(CELLPARMS *parmsPtr, STATE *state, int pOffset, DERIVED *derived, double dt)
{
   PSTATE *pState = (PSTATE *)(((double *)state)+pOffset) ; 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double V = state->Vm; 
   double Cai = state->Cai; 

   double EKs = derived->EKs; 
   double Xs1 = pState->Xs1; 
   double Xs2 = pState->Xs2; 

   double phi =  1 + 0.6*sigm(pow(3.8e-5/Cai,1.4)); 
   derived->I.Ks = cP->GKs *(V-EKs) * phi * Xs1*Xs2;

   double XsMhu = sige(-(V+11.60)/8.932); 
   double Xs1Tau =  817.3 + 1/(2.326e-4*exp((V+48.28)/17.80) + 1.292e-3*exp(-(V+210.0)/230));
   double Xs1TauR =  1/Xs1Tau;
   double Xs2TauR = 1e-2*exp((V-50.0)/20.0) + 1.93e-2*exp(-(V+66.54)/31);
   double dXs1 = (XsMhu-Xs1)*Xs1TauR;  // gate
   double dXs2 = (XsMhu-Xs2)*Xs2TauR;  // gate
   pState->Xs1 += dt*dXs1; 
   pState->Xs2 += dt*dXs2; 
}
// IKs  V,Cai:EKs:GKs:Xs1,Xs2
void OHaraRudy_IKsFuncDebug(CELLPARMS *parmsPtr, STATE *state, int pOffset, DERIVED *derived, double dt)
{
   PSTATE *pState = (PSTATE *)(((double *)state)+pOffset) ; 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double V = state->Vm; 
   double EKs = derived->EKs; 
   double Xs1 = pState->Xs1; 
   double Xs2 = pState->Xs2; 
   double Cai = state->Cai; 

   double phi =  1 + 0.6*sigm(pow(3.8e-5/Cai,1.4)); 
   derived->I.Ks = cP->GKs *(V-EKs) * phi * Xs1*Xs2;

   double XsMhu = sige(-(V+11.60)/8.932); 
   double Xs1Tau =  817.3 + 1/(2.326e-4*exp((V+48.28)/17.80) + 1.292e-3*exp(-(V+210.0)/230));
   double Xs1TauR =  1/Xs1Tau;
   double Xs2TauR = 1e-2*exp((V-50.0)/20.0) + 1.93e-2*exp(-(V+66.54)/31);
   double dXs1 = (XsMhu-Xs1)*Xs1TauR;  // gate
   double dXs2 = (XsMhu-Xs2)*Xs2TauR;  // gate
   pState->Xs1 += dt*dXs1; 
   pState->Xs2 += dt*dXs2; 
}
COMPONENTINFO OHaraRudy_IKsInit()
{
   COMPONENTINFO info;
   info.nVar = nVar; 
   info.varInfo = varInfo;
   info.func = OHaraRudy_IKsFunc;
   info.access = OHaraRudy_IKsAccess;
   return info;
}

