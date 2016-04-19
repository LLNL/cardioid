#include <math.h>
#include "Grandi.h"
#include "Grandi_INa.h"
void Grandi_INaFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt)
{
      VOLTAGE *voltage = (VOLTAGE *)state; 

      PSTATE *pState = (PSTATE *)(state+pOffset) ; 
      PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
      double v = voltage->Vm; 
      double ENa_junc = derived->ENa_junc; 
      double ENa_sl = derived->ENa_sl; 
      double Fsl=1.0-Fjunc;
      
      double m=pState->m; 
      double h=pState->h; 
      double j=pState->j; 

      derived->I.Na_junc=Fjunc*cP->GNa*(v-ENa_junc)*m*m*m*h*j; 
      derived->I.Na_sl=Fsl*cP->GNa*(v-ENa_sl)*m*m*m*h*j; 

      double mss=1.0 / (pow(1.0 + exp( -(56.86 + v) / 9.03 ),2.0));
      double taum=0.1292 * exp(-pow((v+45.79)/15.54,2.0)) + 0.06487 * exp(-pow((v-4.823)/51.12,2.0));
      double ah, bh, aj, bj;
      if (v >= -40.0)
	{
	  ah=0.0;
	  bh=0.77 / (0.13*(1.0 + exp( -(v + 10.66) / 11.1 )));
	  aj=0.0;
	  bj=((0.6 * exp( 0.057 * v)) / (1 + exp( -0.1 * (v + 32) )));  
	}
      else
	{
	  ah=0.057 * exp( -(v + 80.0) / 6.8 );
	  bh=2.7 * exp( 0.079 * v) + 3.1e5 * exp(0.3485 * v);
	  aj=((-2.5428e4*exp(0.2444*v) - 6.948e-6 * exp(-0.04391*v)) * (v + 37.78)) / (1.0 + exp( 0.311 * (v + 79.23) ));
	  bj=(0.02424 * exp( -0.01052 * v )) / (1.0 + exp( -0.1378 * (v + 40.14) ));  
	}

      double tauh=1.0 / (ah + bh);
      double hss=1.0 / (pow(1 + exp( (v + 71.55)/7.43 ),2.0));
      double tauj=1.0 / (aj + bj);
      double jss=1.0 / (pow(1 + exp( (v + 71.55)/7.43 ),2.0));

      double dm=(mss - m) / taum;
      double dh=(hss - h) / tauh;
      double dj=(jss - j) / tauj;

      ENDCODE()

      //pState->m += dt*dm;
      //pState->h += dt*dh;
      //pState->j += dt*dj;

      pState->m=mss-(mss-m)*exp(-dt/taum);
      pState->h=hss-(hss-h)*exp(-dt/tauh);
      pState->j=jss-(jss-j)*exp(-dt/tauj);

}
