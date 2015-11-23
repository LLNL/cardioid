// IKb   V:EK:GKb:
#include <math.h>
#include "OHaraRudy.h"
#include "OHaraRudy_IKb.h"
void OHaraRudy_IKbFunc(CELLPARMS *parmsPtr, double *cell, int pOffset, DERIVED *derived, double dt )
{
   VOLTAGE *voltage = (VOLTAGE *)cell; 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double V = voltage->Vm; 
   double EK = derived->EK;
   double xKb = 1/(1+exp(-(V-14.48)/18.34)) ; 
   derived->I.Kb = cP->GKb*xKb*(V-EK) ;
}
