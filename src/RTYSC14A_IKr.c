// IKr    V::GKr:C3,C2,C1,O,I
#include <math.h>
#include "OHaraRudy.h"
#include "RTYSC14A_IKr.h"
static double a[8]; 
static double b[8]; 
static double TBase = 310.0; // Kelvin

void RTYSC14A_IKrFunc(CELLPARMS *parmsPtr, STATE *state, int pOffset, DERIVED *derived, double dt)
{

   double *S = ((double *)state)+pOffset ; 
   PSTATE *pState = (PSTATE *)S ; 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double Vm = state->Vm; 
   derived->I.Kr = cP->GKr * (Vm-derived->EK)*pState->O; 

   double rate[8]; 
   for (int i=1;i<8;i++) rate[i] = a[i]*exp(b[i]*Vm); 

   double  dSdt[5]; 
   dSdt[0] = rate[1]*S[1] - rate[0]*S[0]; 
   for (int i=1;i<4;i++) 
   {
      dSdt[i] = rate[2*i-2]*S[i-1] + rate[2*i+1]*S[i+1] - (rate[2*i-1]+rate[2*i])*S[i]; 
   } 
   dSdt[4] = rate[6]*S[3] - rate[7]*S[4]; 

   for (int i=0;i<5;i++) S[i] += dt*dSdt[i]; 
}
void RTYSC14A_Rates(double V, double *rate)
{
   rate[0] = T/TBase * exp(24.335 + (T/TBase)*( 0.0112*V-25.914));          //alpha        C3->C2
   rate[1] = T/TBase * exp(13.668 + (TBase/T)*(-0.0603*V-15.707));          //beta         C2->C3
   rate[2] = T/TBase * exp(22.764 + (TBase/T)*(         -25.914));          //alpha_in     C2->C1
   rate[3] = T/TBase * exp(13.193 + (TBase/T)*(         -15.707));          //beta_in      C1->C2
   rate[4] = T/TBase * exp(22.098 + (TBase/T)*( 0.0365*V-25.914));          //alphaalpha   C1->O
   rate[5] = T/TBase * exp( 7.313 + (TBase/T)*(-0.0399*V-15.707));          //betabeta     O->C1
   rate[6] = T/TBase * exp(30.016 + (TBase/T)*( 0.0223*V-30.880))*pow(5.4/Ko,0.4);  // alpha_i O->I
   rate[7] = T/TBase * exp(30.061 + (TBase/T)*(-0.0312*V-33.243));          //beta_i       I->O
}

COMPONENTINFO RTYSC14A_IKrInit()
{
   double rate0[8]; 
   double rate10[8]; 
   RTYSC14A_Rates( 0.0, rate0);
   RTYSC14A_Rates(10.0, rate10);
   for (int i=0;i<8;i++) 
   {
      a[i] = rate0[i];
      b[i] = 0.1*log(rate10[i]/rate0[i]);
   }
   
   COMPONENTINFO info;
   info.nVar = nVar; 
   info.varInfo = varInfo;
   info.func = RTYSC14A_IKrFunc;
   info.access = RTYSC14A_IKrAccess;
   return info;
}
