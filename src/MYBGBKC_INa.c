 
/*
                  BC3  BC2  BC1    BO          10- 11- 12- 13- 14
                   |    |    |     |             
                   C3 - C2 - C1 -  O - IS       0 - 1 - 2 - 3 - 4 
                   |    |    |  /       
                  IC3 -IC2 - IF - IM1-IM2       5 - 6 - 7 - 8 - 9
*/
#include <math.h>
#include "OHaraRudy.h"
#include "MYBGBKC_INa.h"

static double TBase = 310.0; // Kelvin
static double cc; 


void MYBGBKC_Rates(double V, double *rate);

void MYBGBKC_INaFunc(CELLPARMS *parmsPtr, double *cell, int pOffset, DERIVED *derived, double dt)
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
   derived->I.NaL = cP->GNa * (V-derived->ENa)*(pState->O+pState->BO); 
   derived->I.NaFast = 0.0; 

   double rate[50]; 
   MYBGBKC_Rates(V, rate);

   double  dSdt[15]; 
   int nSubSteps=10; 
   for (int l=0;l<nSubSteps;l++)
   {
      for (int i=0;i<15;i++) dSdt[i] = 0.0; 
      for (int k=0;k<3;k++) 
      {
         double *S_k = S + k*5; 
         double *dSdt_k = dSdt + k*5; 
         double *rate_k = rate + k*8; 
         dSdt_k[0] = rate_k[1]*S_k[1] - rate_k[0]*S_k[0]; 
         for (int i=1;i<4;i++) 
         {
            dSdt_k[i] = rate_k[2*i-2]*S_k[i-1] + rate_k[2*i+1]*S_k[i+1] - (rate_k[2*i-1]+rate_k[2*i])*S_k[i]; 
         } 
         dSdt_k[4] = rate_k[6]*S_k[3] - rate_k[7]*S_k[4]; 
      }
      double *rateI = rate+24;
      double *rateB = rate+34;
      double *I = S+5; 
      double *dI = dSdt+5; 
      double *B = S+10; 
      double *dB = dSdt+10; 
      for (int i=0;i<5;i++) 
      {
         dSdt[i] += rateI[2*i]*I[i]+rateB[2*i]*B[i]-(rateI[2*i+1]+rateB[2*i+1])*S[i];
         dI[i]   += rateI[2*i+1]*S[i]   -rateI[2*i]*I[i] ;
         dB[i]   += rateB[2*i+1]*S[i]   -rateB[2*i]*B[i] ;
      }
      dSdt[3] += -rate[44]*S[3] + rate[45]*S[7]; 
      dSdt[7] += +rate[44]*S[3] - rate[45]*S[7]; 
      if (l==0) 
      {
         if (derived->dState != 0) 
         {
               double  *dState = derived->dState+pOffset;
               for (int i=0;i<15;i++) dState[i] = dSdt[i]; 
         }
      }
      for (int i=0;i<15;i++) S[i] += dt*dSdt[i]/nSubSteps; 
   }
      


   /*
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
      dSdt[3] += cP->rC*S[8]   ; dSdt[8] -= cP->rC*S[8]; 
      dSdt[4] += cP->rC*S[9]   ; dSdt[9] -= cP->rC*S[9]; 
    */


}
void MYBGBKC_Rates(double V, double *rate)
{

   double  aa[] = {-1/8.2471,-1/12.0,-1/12.542,-1/15.0,-1/17.0,-1/20.3,-1/43.749,-1/150};
   double  ex[8]; 
   for (int i=0;i<46;i++) rate[i] = 0; 

   double alpha11 = 8.5539/(7.4392e-2*exp(-V/17.0)+2.0373e-1*exp(-V/150));   //alpha11 C3->C2, IC3->IC2
   double alpha12 = 8.5539/(7.4392e-2*exp(-V/15.0)+2.0373e-1*exp(-V/150));   //alpha12 C2->C1, IC2->IF
   double alpha13 = 8.5539/(7.4392e-2*exp(-V/12.0)+2.0373e-1*exp(-V/150));   //alpha13 C1 ->O
   double beta11  = 7.5215e-2*exp(-(V   )/20.3);                             //beta11  C2->C3  IC2->IC3
   double beta12  = 2.7574   *exp(-(V- 5)/20.3);                             //beta12  C1->C2, IF->IC2
   double beta13  = 4.7755e-1*exp(-(V-10)/20.3);                             //beta13  O->C1
   double alpha3  = 5.1458e-6*exp(-V/8.2471);                                //alpha3  IC3->C3, IC2->C2, IF->C1
   double beta3   = 6.1205   *exp( V/13.542);      //12.542 is paper value   //beta3   C3->IC3, C2->IC2, C1->IF
   double alpha2  = 13.370   *exp( V/43.749);                                //alpha2  O->IF
   double beta2   = (alpha13*alpha2*alpha3)/(beta13*beta3);                  //beta2   IF->O
   double alphax  = 3.4229e-2*alpha2;                                        //alphax  O->IS
   double betax   = 1.7898e-2*alpha3;                                        //betax   IS->O


   double mu1WT   = 2.0462e-7;                                               //mu1WT   C3, C2, C1, O -> BC3, BC2, BC1, BO  
   double mu1HF   = 2.7252e-7;                                               //mu1HF   C3, C2, C1, O -> BC3, BC2, BC1, BO
   double mu2WT   = 8.9731e-4;                                               //mu2WT   BC3, BC2, BC1, BO -> C3, C2, C1, O
   double mu2HF   = 1.9701e-4;                                               //mu2HF   BC3, BC2, BC1, BO -> C3, C2, C1, O


   int k=0; 
   rate[k++]  =  alpha11;
   rate[k++]  =  beta11;
   rate[k++]  =  alpha12;
   rate[k++]  =  beta12;
   rate[k++]  =  alpha13;
   rate[k++]  =  beta13;
   rate[k++]  =  alphax ;
   rate[k++]  =  betax ;

   rate[k++] =  alpha11;
   rate[k++] =  beta11;
   rate[k++] =  alpha12;
   rate[k++] =  beta12;
   rate[k++] =  0.0; 
   rate[k++] =  0.0; 
   rate[k++] =  0.0; 
   rate[k++] =  0.0; 

   rate[k++] =  alpha11;
   rate[k++] =  beta11;
   rate[k++] =  alpha12;
   rate[k++] =  beta12;
   rate[k++] =  alpha13; 
   rate[k++] =  beta13;
   rate[k++] =  0.0; 
   rate[k++] =  0.0; 


   rate[k++] =  alpha3;
   rate[k++] =  beta3;
   rate[k++] =  alpha3;
   rate[k++] =  beta3;
   rate[k++] =  alpha3; 
   rate[k++] =  beta3; 
   rate[k++] =  0.0;   
   rate[k++] =  0.0;   
   rate[k++] =  0.0; 
   rate[k++] =  0.0; 

   rate[k++] = mu2WT; 
   rate[k++] = mu1WT; 
   rate[k++] = mu2WT; 
   rate[k++] = mu1WT; 
   rate[k++] = mu2WT; 
   rate[k++] = mu1WT; 
   rate[k++] = mu2WT; 
   rate[k++] = mu1WT; 
   rate[k++] = 0.0; 
   rate[k++] = 0.0; 

   rate[k++]  =  alpha2;
   rate[k++]  =  beta2;


}
void MYBGBKC_INaConstants()
{
   cc = sqrt(Ko/5.4); 
}
