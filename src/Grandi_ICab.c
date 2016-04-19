#include <math.h>
#include "Grandi.h"
#include "Grandi_ICab.h"

void Grandi_ICabFunc(CELLPARMS *parmsPtr, double *cell, int pOffset, DERIVED *derived,  double dt)
{
   VOLTAGE *voltage = (VOLTAGE *)cell; 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double v = voltage->Vm; 
   double ECa_junc = derived->ECa_junc; 
   double ECa_sl = derived->ECa_sl; 
   double Fsl=1.0-Fjunc;

   derived->I.Cab_junc=Fjunc*cP->GCab*(v-ECa_junc);
   derived->I.Cab_sl=Fsl*cP->GCab*(v-ECa_sl);
}
