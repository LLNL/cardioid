// INaFast V:ENa,phiCaMK:GNaFast:m,hFast,hSlow,j,hCaMKSlow,jCaMK
#include <math.h>
#include "OHaraRudy.h"
#include "OHaraRudy_INaFast.h"
void OHaraRudy_INaFastFunc(CELLPARMS *parmsPtr, STATE *state, int pOffset, DERIVED *derived, double dt)
{
//#define hSlowTauC   0.009764   //paper
#define hSlowTauC   0.009794   //code

   PSTATE *pState = (PSTATE *)(((double *)state)+pOffset) ; 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double V = state->Vm; 
   double ENa = derived->ENa; 
   double phiCaMK=derived->phiCaMK;

      //Gates needed to calculate INaFast; 
      double m=pState->m; 
      double hFast=pState->hFast; 
      double hSlow=pState->hSlow; 
      double j=pState->j; 
      double hCaMKSlow=pState->hCaMKSlow;
      double jCaMK = pState->jCaMK;   

      double AhFast=0.99;
      double AhSlow = 1-AhFast; 
      double h     = AhFast*hFast     + AhSlow*hSlow; 
      double AhCaMKFast = AhFast; 
      double AhCaMKSlow = AhSlow; 
      double hCaMKFast = hFast; 
      double hCaMK = AhCaMKFast*hCaMKFast + AhCaMKSlow*hCaMKSlow; 

      derived->I.NaFast = cP->GNaFast*(V-ENa)*m*m*m*((1-phiCaMK)*h*j+phiCaMK*hCaMK*jCaMK); 

      double mMhu =  sige(-(V + 39.57)/9.871);  //OHR orginal  mMhu
//    double mMhu = SQ(sige(-(V + 56.86)/9.030)); //TT06  mMhu
      double mTauR = 6.765*exp((V + 11.64)/34.77)+8.552*exp(-(V + 77.42)/5.955); 
      double  dm = (mMhu-m)*mTauR;  // gate

      double hMhu = sige((V + 82.90)/6.086); 
      double hFastTauR = 1.432e-5 * exp(-(V + 1.196)/6.285) + 6.149 * exp((V + 0.5096)/20.27);
      double hSlowTauR = hSlowTauC * exp (-(V + 17.95)/28.05) + 0.3343 * exp((V+5.730)/56.66);
      double dhFast = (hMhu - hFast)* hFastTauR;  // gate
      double dhSlow = (hMhu - hSlow)* hSlowTauR;  // gate

      double jMhu = hMhu; 
      double jTau  = 2.038 + 1.0/(0.02136*exp(-(V + 100.6)/8.281) + 0.3052*exp((V + 0.9941)/38.45));
      double jTauR = 1.0/jTau; 
      double dj = (jMhu-j)*jTauR;  // gate

      double jCaMKMhu = jMhu; 
      double jCaMKTauR = jTauR/1.46; 
      double djCaMK = (jCaMKMhu-jCaMK)*jCaMKTauR;  // gate

      double hCaMKMhu = sige((V+89.1)/6.086);
      double hCaMKSlowTauR = hSlowTauR/3.0; 
      double dhCaMKSlow = (hCaMKMhu-hCaMKSlow)*hCaMKSlowTauR;  // gate
      pState->m += dt*dm;    //Forward Euler for original OR mGate 
//    double tauRdt = 1-exp(-dt*mTauR); //Rush Larsen needed with TT06 mMhu; 
//    pState->m += (mMhu-m)*tauRdt;     // Rush Larsen needed with TT06 mMhu; 
      pState->hFast += dt*dhFast; 
      pState->hSlow += dt*dhSlow; 
      pState->j     += dt*dj; 
      pState->hCaMKSlow += dt*dhCaMKSlow; 
      pState->jCaMK += dt*djCaMK;
}
