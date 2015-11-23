// INaCai   V,Cai,Nai::GNaCai:
// INaCass  V,Cass,Nass::GNaCass:
#include "OHaraRudy.h"
#include "OHaraRudy_INaCaCommon.h"
#include "OHaraRudy_INaCai.h"
void OHaraRudy_INaCaiFunc(CELLPARMS *parmsPtr, double *cell, int pOffset, DERIVED *derived, double dt)
{
   VOLTAGE *voltage = (VOLTAGE *)cell; 
   CONCENTRATIONS   *concentrations = (CONCENTRATIONS*) (cell + CONCENTRATIONS_OFFSET); 
   double Cai = concentrations->Cai; 
   double Nai = concentrations->Nai; 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double V = voltage->Vm; 
   derived->I.NaCai =  cP->GNaCai *INaCaYFunc(V,Cai, Nai   );
}
