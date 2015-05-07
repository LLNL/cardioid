#include "OHaraRudy.h"
#include <math.h>

static double RTF ;
static double RTFlogNao; 
static double RTFlogKo; 
static double RTFlogKoNao; 

void reversalPotentialsInit()
{
   RTF = (R*T)/F;
   RTFlogNao= RTF*log(Nao); 
   RTFlogKo= RTF*log(Ko); 
   RTFlogKoNao=RTF*log(Ko+PRNaK*Nao);
}
void reversalPotentials(double Nai, double Ki,DERIVED *derived)
{
   derived->ENa = RTFlogNao-RTF*log(Nai);
   derived->EK  = RTFlogKo -RTF*log(Ki);
   derived->EKs = RTFlogKoNao-RTF*log((Ki+PRNaK*Nai));
}
