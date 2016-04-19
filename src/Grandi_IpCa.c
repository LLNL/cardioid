#include <math.h>
#include "Grandi.h"
#include "Grandi_IpCa.h"

void Grandi_IpCaFunc(CELLPARMS *parmsPtr, double *cell, int pOffset, DERIVED *derived, double dt )
{
   VOLTAGE *voltage = (VOLTAGE *)cell; 
   CONCENTRATIONS   *concentrations = (CONCENTRATIONS*) (cell + CONCENTRATIONS_OFFSET); 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double Caj = concentrations->Caj; 
   double Casl = concentrations->Casl;

   double Qpow=(T-310.0)/10.0;
   double Q10SLCaP=2.35;
   double  KmPCa=0.5e-3;
   double Fsl=1.0-Fjunc;

   derived->I.pCa_junc=Fjunc*pow(Q10SLCaP,Qpow)*cP->IbarSLCaP*pow(Caj,1.6)/(pow(KmPCa,1.6)+pow(Caj,1.6));
   derived->I.pCa_sl=Fsl*pow(Q10SLCaP,Qpow)*cP->IbarSLCaP*pow(Casl,1.6)/(pow(KmPCa,1.6)+pow(Casl,1.6));
}
