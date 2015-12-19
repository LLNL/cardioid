// IKr    V::GKr:C3,C2,C1,O,I
#include <math.h>
#include "OHaraRudy.h"
#include "RTYSC14A_IKr.h"
static double a[8]; 
static double b[8]; 
static double TBase = 310.0; // Kelvin
static double cc; 

void RTYSC14A_IKrFunc(CELLPARMS *parmsPtr, double *cell, int pOffset, DERIVED *derived, double dt)
{

   VOLTAGE *voltage = (VOLTAGE *)cell; 
   double *S = cell+pOffset ; 
   PSTATE *pState = (PSTATE *)S ; 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
/*
   cP->rC = 0.00;  
   cP->kC = 0.00;  
   cP->rO = 16.227;
   cP->kO =  0.511;
   cP->rI =  0.2318;
   cP->kI =  0.511;
// cP->D = 48e-6; 
   cP->D = 0.0; 
   cP->rO = 0.0;    
   cP->kO =  0.0;    
   cP->rI =  0.0;    
   cP->kI =  0.0;    
   cP->D = 0.0;    
*/
   double V = voltage->Vm; 
   derived->I.Kr = cP->GKr * cc* (V-derived->EK)*pState->O; 

   double rate[8]; 
   for (int i=0;i<8;i++) rate[i] = a[i]*exp(b[i]*V); 

   double  dSdt[10]; 
   dSdt[0] = rate[1]*S[1] - rate[0]*S[0]; 
   for (int i=1;i<4;i++) 
   {
      dSdt[i] = rate[2*i-2]*S[i-1] + rate[2*i+1]*S[i+1] - (rate[2*i-1]+rate[2*i])*S[i]; 
   } 
   dSdt[4] = rate[6]*S[3] - rate[7]*S[4]; 

   dSdt[0+5] = rate[1]*S[1+5] - rate[0]*S[0+5]; 
   for (int i=1;i<4;i++) 
   {
      dSdt[i+5] = rate[2*i-2]*S[i-1+5] + rate[2*i+1]*S[i+1+5] - (rate[2*i-1]+rate[2*i])*S[i+5]; 
   } 
   dSdt[4+5] = rate[6]*S[3+5] - rate[7]*S[4+5]; 

   double D = cP->D; 
   dSdt[0] -= cP->kC*D*S[0] ; dSdt[5] += cP->kC*D*S[0]; 
   dSdt[1] -= cP->kC*D*S[1] ; dSdt[6] += cP->kC*D*S[1]; 
   dSdt[2] -= cP->kC*D*S[2] ; dSdt[7] += cP->kC*D*S[2]; 
   dSdt[3] -= cP->kO*D*S[3] ; dSdt[8] += cP->kO*D*S[3]; 
   dSdt[4] -= cP->kI*D*S[4] ; dSdt[9] += cP->kI*D*S[4]; 
   dSdt[0] += cP->rC*S[5]   ; dSdt[5] -= cP->rC*S[5]; 
   dSdt[1] += cP->rC*S[6]   ; dSdt[6] -= cP->rC*S[6]; 
   dSdt[2] += cP->rC*S[7]   ; dSdt[7] -= cP->rC*S[7]; 
   dSdt[3] += cP->rO*S[8]   ; dSdt[8] -= cP->rO*S[8]; 
   dSdt[4] += cP->rI*S[9]   ; dSdt[9] -= cP->rI*S[9]; 
   if (derived->dState != 0) 
   {
      double  *dState = derived->dState+pOffset;
      for (int i=0;i<10;i++) dState[i] = dSdt[i]; 
   }
   for (int i=0;i<10;i++) S[i] += dt*dSdt[i]; 

}
void RTYSC14A_Rates(double V, double *rate)
{
   rate[0] = T/TBase * exp(24.335 + (T/TBase)*( 0.0112*V-25.914));                  //alpha        C3->C2
   rate[1] = T/TBase * exp(13.688 + (TBase/T)*(-0.0603*V-15.707));                  //beta         C2->C3
   rate[2] = T/TBase * exp(22.746 + (TBase/T)*(         -25.914));                  //alpha_in     C2->C1
   rate[3] = T/TBase * exp(13.193 + (TBase/T)*(         -15.707));                  //beta_in      C1->C2
   rate[4] = T/TBase * exp(22.098 + (TBase/T)*( 0.0365*V-25.914));                  //alphaalpha   C1->O
   rate[5] = T/TBase * exp( 7.313 + (TBase/T)*(-0.0399*V-15.707));                  //betabeta     O->C1
   rate[6] = T/TBase * exp(30.016 + (TBase/T)*( 0.0223*V-30.888))*pow(5.4/Ko,0.4);  // alpha_i     O->I  // paper used 30.880
   rate[7] = T/TBase * exp(30.061 + (TBase/T)*(-0.0312*V-33.243));                  //beta_i       I->O
}

void RTYSC14A_IKrConstants()
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
   cc = sqrt(Ko/5.4); 
}
