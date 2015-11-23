// INaCai   V,Cai,Nai::GNaCai:
// INaCass  V,Cass,Nass::GNaCass:
#include "OHaraRudy.h"
#include "OHaraRudy_INaCaCommon.h"
#include "OHaraRudy_INaCass.h"


void OHaraRudy_INaCassFunc(CELLPARMS *parmsPtr, double *cell, int pOffset, DERIVED *derived, double dt)
{
   VOLTAGE *voltage = (VOLTAGE *)cell; 
   CONCENTRATIONS   *concentrations = (CONCENTRATIONS*) (cell + CONCENTRATIONS_OFFSET); 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double V = voltage->Vm; 
   double Cass = concentrations->Cass; 
   double Nass = concentrations->Nass; 
   derived->I.NaCass = cP->GNaCass*INaCaYFunc(V,Cass, Nass  );
}
