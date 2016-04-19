#include "Grandi.h"
#include <math.h>

static double RTF;
static double FRT;

void Grandi_RevpotsInit()
{
   RTF = (R*T)/F;
}
void Grandi_Revpots(double Naj, double Nasl, double Nai, double Ki, double Caj, double Casl, DERIVED *derived)
{
   derived->ENa_junc=RTF*log(Nao/Naj);     // [mV]
   derived->ENa_sl=RTF*log(Nao/Nasl);       // [mV]
   derived->EK=RTF*log(Ko/Ki);	        // [mV]
   derived->EKs=RTF*log((Ko+pNaK*Nao)/(Ki+pNaK*Nai));
   derived->ECa_junc=(RTF/2.0)*log(Cao/Caj);   // [mV]
   derived->ECa_sl=(RTF/2.0)*log(Cao/Casl);     // [mV]
   derived->ECl=RTF*log(Cli/Clo);            // [mV]
}
