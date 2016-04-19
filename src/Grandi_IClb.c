#include <math.h>
#include "Grandi.h"
#include "Grandi_IClb.h"

void Grandi_IClbFunc(CELLPARMS *parmsPtr, double *cell, int pOffset, DERIVED *derived, double dt )
{
#define KdClb 100.0e-3;

   VOLTAGE *voltage = (VOLTAGE *)cell; 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double v = voltage->Vm; 
   double ECl = derived->ECl;

   derived->I.Clb=cP->GClb*(v-ECl);
}
