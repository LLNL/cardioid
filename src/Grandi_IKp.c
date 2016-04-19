#include <math.h>
#include "Grandi.h"
#include "Grandi_IKp.h"
void Grandi_IKpFunc(CELLPARMS *parmsPtr, double *cell, int pOffset, DERIVED *derived, double dt )
{
   VOLTAGE *voltage = (VOLTAGE *)cell; 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double v = voltage->Vm; 
   double EK = derived->EK;
   double Fsl=1.0-Fjunc;

   double xkp = 1.0/(1.0+exp(7.488-v/5.98));

   derived->I.Kp_junc=Fjunc*cP->GKp*xkp*(v-EK);
   derived->I.Kp_sl=Fsl*cP->GKp*xkp*(v-EK);
}
