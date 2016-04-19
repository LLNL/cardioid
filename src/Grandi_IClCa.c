#include <math.h>
#include "Grandi.h"
#include "Grandi_IClCa.h"

void Grandi_IClCaFunc(CELLPARMS *parmsPtr, double *cell, int pOffset, DERIVED *derived, double dt )
{

   VOLTAGE *voltage = (VOLTAGE *)cell; 
   CONCENTRATIONS   *concentrations = (CONCENTRATIONS*) (cell + CONCENTRATIONS_OFFSET); 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double v = voltage->Vm; 
   double ECl = derived->ECl;
   double Caj=concentrations->Caj;
   double Casl=concentrations->Casl;

   double KdClCa=100.0e-3;
   double Fsl=1.0-Fjunc;

   derived->I.ClCa_junc=Fjunc*cP->GClCa/(1.0+KdClCa/Caj)*(v-ECl);
   derived->I.ClCa_sl=Fsl*cP->GClCa/(1.0+KdClCa/Casl)*(v-ECl);
}
