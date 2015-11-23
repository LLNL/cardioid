// Ito V:EK,phiCaMK,:Gto:a,iFast,iSlow,aCaMK,iCaMKFast,iCaMKSlow
#include <math.h>
#include "OHaraRudy.h"
#include "OHaraRudy_Ito.h"

void OHaraRudy_ItoFunc(CELLPARMS *parmsPtr, double *cell, int pOffset, DERIVED *derived, double dt )
{
#define iSlowTauC   1.780e-8 // code
#define aTauC0     18.4099   // code
#define aTauC1     29.3814   // code

   //#define iSlowTauC 1.7808e-8 // paper
   //#define aTauC0   18.41      // paper 
   //#define aTauC1   29.38      // paper 

   VOLTAGE *voltage = (VOLTAGE *)cell; 
   PSTATE *pState = (PSTATE *)(((double *)cell)+pOffset) ; 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double V = voltage->Vm; 
   double EK = derived->EK; 
   double phiCaMK   = derived->phiCaMK; 

   double a     =pState->a;
   double iFast = pState->iFast; 
   double iSlow = pState->iSlow; 
   double aCaMK = pState->aCaMK;
   double iCaMKFast = pState->iCaMKFast; 
   double iCaMKSlow = pState->iCaMKSlow; 

   double AiFast = sige((V-213.6)/151.2); 
   double AiSlow = 1 - AiFast; 
   double AiCaMKFast = AiFast; 
   double AiCaMKSlow = AiSlow; 

   double i = AiFast * iFast + AiSlow * iSlow; 
   double iCaMK = AiCaMKFast * iCaMKFast + AiCaMKSlow * iCaMKSlow; 

   derived->I.to=cP->Gto*(V-EK) * ((1-phiCaMK)*a*i + phiCaMK *aCaMK*iCaMK); 

   double aMhu = sige(-(V-14.34)/14.82); 
 //double aTauR = ((1./1.2089)*sige(-(V-18.4099)/29.3814) + 3.5*sige((V+100)/29.3814))/1.0515;
   double aTauR = ((1./1.2089)*sige(-(V-aTauC0)/29.3814) + 3.5*sige((V+100)/aTauC1))/1.0515;
   double da = (aMhu-a)*aTauR;  // gate

   double iMhu = sige((V+43.94)/5.711); 
   double delta = 1.0-cP->aDelta*sige((V+70)/5.0); 
   double iFastTau = 4.562 + 1/(0.3933*exp(-(V+100)/100) + 0.08004*exp((V+50)/16.59));
   iFastTau = iFastTau*delta; 
   //double iSlowTau = 23.62 + 1/(0.001416*exp(-(V+96.52)/59.05) + 1.7808e-8*exp((V+114.1)/8.079));
   double iSlowTau = 23.62 + 1/(0.001416*exp(-(V+96.52)/59.05) + iSlowTauC*exp((V+114.1)/8.079));
   iSlowTau = iSlowTau*delta; 
   double iFastTauR = 1/iFastTau; 
   double iSlowTauR = 1/iSlowTau; 

   double diSlow = (iMhu-iSlow)*iSlowTauR;  // gate
   double diFast = (iMhu-iFast)*iFastTauR;  // gate

   double aCaMKMhu = sige(-(V-24.34)/14.82); 
   double aCaMKTauR = aTauR; 
   double daCaMK = (aCaMKMhu-aCaMK)*aCaMKTauR;  // gate

   double iCaMKMhu = iMhu; 
   double deltaCaMKdevelop = 1.354 + 1e-4/(exp((V-167.4)/15.89)+exp(-(V-12.23)/0.2154));
   double deltaCaMKrecover = 1.0 - 0.5*sige((V+70)/20);
   double deltaCaMKR =  1/(deltaCaMKdevelop*deltaCaMKrecover); 
   double iCaMKFastTauR = iFastTauR * deltaCaMKR;
   double iCaMKSlowTauR = iSlowTauR * deltaCaMKR;
   double diCaMKSlow = (iCaMKMhu-iCaMKSlow)*iCaMKSlowTauR;  // gate
   double diCaMKFast = (iCaMKMhu-iCaMKFast)*iCaMKFastTauR;  // gate
   ENDCODE()
   pState->a += dt*da; 
   pState->iSlow += dt*diSlow; 
   pState->iFast += dt*diFast; 
   pState->aCaMK += dt*daCaMK; 
   pState->iCaMKSlow += dt*diCaMKSlow; 
   pState->iCaMKFast += dt*diCaMKFast; 
}
