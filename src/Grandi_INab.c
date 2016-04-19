#include <math.h>
#include "Grandi.h"
#include "Grandi_INab.h"

void Grandi_INabFunc(CELLPARMS *parmsPtr, double *cell, int pOffset, DERIVED *derived,  double dt)
{
   VOLTAGE *voltage = (VOLTAGE *)cell; 
   CONCENTRATIONS   *concentrations = (CONCENTRATIONS*) (cell + CONCENTRATIONS_OFFSET); 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double v = voltage->Vm; 
   double ENa_junc = derived->ENa_junc; 
   double ENa_sl = derived->ENa_sl; 
   double Fsl=1.0-Fjunc;

   derived->I.Nab_junc=Fjunc*cP->GNab*(v-ENa_junc);
   derived->I.Nab_sl=Fsl*cP->GNab*(v-ENa_sl);
}
