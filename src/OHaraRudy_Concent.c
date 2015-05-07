#include <math.h>
#include "OHaraRudy.h"
#include "OHaraRudy_Concent.h"
static double TRPN =  0.07 ;  // mM
static double BSR  =  0.047;  // mM
static double BSL  =  1.124;  // mM
static double CSQN = 10.0  ;  // mM
static double KmCMDN = 0.00238; // mM
static double KmTRPN = 0.0005 ; // mM
static double KmBSR  = 0.00087; // mM
static double KmBSL  = 0.0087 ; // mM
static double KmCSQN = 0.8    ; // mM


//Cell Geometry Parameters 
static   double AFV     ;//= Acap/(F*Vmyo);
static   double AFVss   ;//= Acap/(F*Vss);
static   double VssVmyo ;// = Vss/Vmyo; 
static   double VnsrVmyo;//= Vnsr/Vmyo; 
static   double VjsrVss ;//= Vjsr/Vss; 
static   double VjsrVnsr;//= Vjsr/Vnsr; 

void OHaraRudyCellular()
{
   double L = 0.01; // cm 
   double r = 0.0011; // cm
   double RCG = 2.0; 
   double Vcell ; // uL   Vcell = pi*r^2*L  Volume of cell; 
   double Ageo  ;// cm^2  Ageo =  2*pi*(r^2 + r*L)  Surface Area of cell 
   double Acap  ;// cm^2  Acap = RCG*Ageo 
   double Vmyo  ;// uL   Vmyo = 0.6800*Vcell 
   double Vnsr  ;// uL   Vnsr = 0.0552*Vcell 
   double Vjsr  ;// uL   Vjsr = 0.0048*Vcell 
   double Vss   ;// uL   Vss  = 0.0200*Vcell 
   double pi = 3.14; 
   Vcell = pi *r*r*L * 1e3; 
   Ageo  = 2*pi*(r*r+r*L);
   Acap  = RCG*Ageo; 
   Vmyo = 0.68*Vcell; 
   Vnsr = 0.0552*Vcell; 
   Vjsr = 0.0048*Vcell; 
   Vss = 0.02*Vcell; 
//
   AFV = Acap/(F*Vmyo);
   AFVss = Acap/(F*Vss);
   VssVmyo  = Vss/Vmyo; 
   VnsrVmyo = Vnsr/Vmyo; 
   VjsrVss  = Vjsr/Vss; 
   VjsrVnsr  = Vjsr/Vnsr; 
}
void OHaraRudy_ConcentFunc(CELLPARMS *parmsPtr, STATE *state, int pOffset, DERIVED *derived, double dt)
{

   PSTATE *pState = (PSTATE *)(((double *)state)+pOffset) ; 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double V    = pState->Vm; 
   double Nass = pState->Nass; 
   double Nai  = pState->Nai; 
   double Cass = pState->Cass; 
   double Cajsr= pState->Cajsr; 
   double Cansr= pState->Cansr; 
   double Cai  = pState->Cai; 
   double Kss  = pState->Kss; 
   double Ki   = pState->Ki; 

   CURRENTS *I = &(derived->I); 
// Fluxes 
   FLUXES J =derived->J; 

   double INa = I->NaFast + I->NaL + I->Nab; 
   double IK  = I->to + I->Kr + I->Ks + I->K1 + I->Kb;
   double IpCab = I->pCa + I->Cab; 

   double bCai   = 1/(1 + (cP->CMDN/KmCMDN)*sigm2(Cai/KmCMDN)+(TRPN/KmTRPN)*sigm2(Cai/KmTRPN));
   double bCass  = 1/(1 + (BSR /KmBSR )*sigm2(Cass/KmBSR)+(BSL/KmBSL)*sigm2(Cass/KmBSL));
   double bCajsr = 1/(1+(CSQN/KmCSQN)*sigm2(Cajsr/KmCSQN)); 

   double dNai = -(INa + 3*I->NaCai + 3*I->NaK )*AFV            + J.diffNa*VssVmyo;
   double dNass= -(I->CaNa + 3*I->NaCass       )*AFVss          - J.diffNa; 
   double dKi =  -(IK - 2*I->NaK + I->stimulus       )*AFV            + J.diffK*VssVmyo;
   double dKss = -(I->CaK                      )*AFVss          - J.diffK; 
   double dCai= (-(IpCab  - 2*I->NaCai         )*AFV  *0.5      - J.up*VnsrVmyo + J.diffCa*VssVmyo)*bCai; 
   double dCass=(-(I->CaL - 2*I->NaCass        )*AFVss*0.5      + J.rel*VjsrVss - J.diffCa        )*bCass; 
   double dCansr =                                                J.up-J.tr*VjsrVnsr; 
   double dCajsr =                                               (J.tr-J.rel                      )*bCajsr; 
   double dVm = -(INa + IK + I->CaL + I->CaNa + I->CaK + I->NaCai + I->NaCass + I->NaK + IpCab )/Cm;
   derived->dVm = dVm; 

   pState->Nai  *= exp(dt*dNai/Nai); 
   pState->Nass *= exp(dt*dNass/Nass); 
   pState->Ki   *= exp(dt*dKi/Ki); 
   pState->Kss  *= exp(dt*dKss/Kss); 
   pState->Cai  *= exp(dt*dCai/Cai); 
   pState->Cass *= exp(dt*dCass/Cass); 
   pState->Cansr*= exp(dt*dCansr/state->Cansr); 
   pState->Cajsr*= exp(dt*dCajsr/state->Cajsr); 

}
COMPONENTINFO OHaraRudy_ConcentInit()
{
   OHaraRudyCellular();
   COMPONENTINFO info;
   info.nVar = nVar; 
   info.varInfo = varInfo;
   info.func = OHaraRudy_ConcentFunc;
   info.access = OHaraRudy_ConcentAccess;
   return info;
}
