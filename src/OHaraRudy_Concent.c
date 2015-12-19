#include <math.h>
#include <stdlib.h>
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
static   double Acap;  ;
static   double AFV     ;//= Acap/(F*Vmyo);
static   double AFVss   ;//= Acap/(F*Vss);
static   double VssVmyo ;// = Vss/Vmyo; 
static   double VnsrVmyo;//= Vnsr/Vmyo; 
static   double VjsrVss ;//= Vjsr/Vss; 
static   double VjsrVnsr;//= Vjsr/Vnsr; 
static   double Vmyo  ;// uL   Vmyo = 0.6800*Vcell 
static   double Vnsr  ;// uL   Vnsr = 0.0552*Vcell 
static   double Vjsr  ;// uL   Vjsr = 0.0048*Vcell 
static   double Vss   ;// uL   Vss  = 0.0200*Vcell 

void OHaraRudy_ConcentConstants()
{
   double pi = 3.14; 
   double L = 0.01; // cm 
   double r = 0.0011; // cm
   double RCG = 2.0; 
   double Vcell ; // uL   Vcell = pi*r^2*L  Volume of cell; 
   double Ageo  ;// cm^2  Ageo =  2*pi*(r^2 + r*L)  Surface Area of cell 
   Acap  ;// cm^2  Acap = RCG*Ageo 
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
void OHaraRudy_ConcentFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt)
{
   VOLTAGE *voltage = (VOLTAGE *)state; 
   PSTATE *pState = (PSTATE *)(state+pOffset) ; 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   //double V = voltage->Vm; 
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
   I->stimulus = voltage->iStim; 

   double INa = I->NaFast + I->NaL + I->Nab; 
   double IK  = I->to + I->Kr + I->Ks + I->K1 + I->Kb;
   double IpCab = I->pCa + I->Cab; 

   double bCai   = 1/(1 + (cP->CMDN/KmCMDN)*sigm2(Cai/KmCMDN)+(TRPN/KmTRPN)*sigm2(Cai/KmTRPN));
   double bCass  = 1/(1 + (BSR /KmBSR )*sigm2(Cass/KmBSR)+(BSL/KmBSL)*sigm2(Cass/KmBSL));
   double bCajsr = 1/(1+(CSQN/KmCSQN)*sigm2(Cajsr/KmCSQN)); 


   double dNai = -(INa + 3*I->NaCai + 3*I->NaK )*AFV            + J.diffNa*VssVmyo;

   double dNass= -(I->CaNa + 3*I->NaCass       )*AFVss          - J.diffNa; 
   double dKi =  -(IK - 2*I->NaK + I->stimulus )*AFV            + J.diffK*VssVmyo;                 //             dVm = dVmR + dVmD + dVmS

   double dKss = -(I->CaK                      )*AFVss          - J.diffK; 
   double dCai= (-(IpCab  - 2*I->NaCai         )*AFV  *0.5      - J.up*VnsrVmyo + J.diffCa*VssVmyo)*bCai; 
   double dCass=(-(I->CaL - 2*I->NaCass        )*AFVss*0.5      + J.rel*VjsrVss - J.diffCa        )*bCass; 
   double dCansr =                                                J.up-J.tr*VjsrVnsr; 
   double dCajsr =                                               (J.tr-J.rel                      )*bCajsr; 
   double dVm = -(INa + IK + I->CaL + I->CaNa + I->CaK + I->NaCai + I->NaCass + I->NaK + IpCab )/Cm;
//printf("%24.15e %24.15e\n",dVm,Cm); 
//printf("%24.15e\n%24.15e\n%24.15e\n%24.15e\n%24.15e\n%24.15e\n%24.15e\n%24.15e\n%24.15e\n%24.15e\n%24.15e\n%24.15e\n%24.15e\n%24.15e\n%24.15e\n",
//        I->NaFast,I->NaL,I->to,I->CaL,I->CaNa,I->CaK,I->Kr,I->Ks,I->K1,I->NaCai+I->NaCass,I->NaK,I->Nab,I->Kb,I->pCa,I->Cab);

   voltage->dVm = dVm; 
   derived->dState[0] = dVm; 
//#ifdef SA
   static int loop = 0; 
   //double CaiBound  = cP->CMDN * (1-KmCMDN/(KmCMDN + Cai))  + TRPN*(1-KmTRPN/(KmTRPN + Cai));    //mM
   //double CassBound =      BSR * (1-KmBSR /(KmBSR + Cass)) +  BSL*(1-KmBSL /(KmBSL  + Cass));    //mM
   //double CajsrBound=     CSQN * (1-KmCSQN/(KmCSQN + Cajsr));   //mM
   //double CaTotal = (Cai + CaiBound)*Vmyo + (Cass + CassBound)*Vss + (Cajsr + CajsrBound)*Vjsr + Cansr*Vnsr;  //milli-mole
   //double NaTotal = Nai*Vmyo + Nass*Vss;     //10^-9 *mole
   //double KTotal = Ki*Vmyo  + Kss*Vss;      // 10^-9 *mole
   //double qTotal = (2*CaTotal+NaTotal+KTotal)*F;  //10^-9 C
   //double CTotal = Cm*Acap; // uF
   //printf("qTotal = %e (10^-3 uC) \n",qTotal); 
   //printf("CTotal = %e (uF) \n",CTotal); 
   //static double qo;
   //double V1 = (qTotal-qo)/CTotal; 
   //if (loop == 0) qo = (qTotal - V*CTotal);
   //if ((loop-1) % 100000 ==0) printf("%8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",dNai/Nai,dNass/Nass,dKi/Ki,dKss/Kss,dCai/Cai,dCass/Cass,dCansr/Cansr,dCajsr/Cajsr,I->NaFast,I->NaL,I->Nab,I->NaCai,I->NaK,I->NaCass,I->CaNa); 

   loop++; 
//#endif

   ENDCODE()
   pState->Nai  *= exp(dt*dNai/Nai); 
   pState->Nass *= exp(dt*dNass/Nass); 
   pState->Ki   *= exp(dt*dKi/Ki); 
   pState->Kss  *= exp(dt*dKss/Kss); 
   pState->Cai  *= exp(dt*dCai/Cai); 
   pState->Cass *= exp(dt*dCass/Cass); 
   pState->Cansr*= exp(dt*dCansr/Cansr); 
   pState->Cajsr*= exp(dt*dCajsr/Cajsr); 

}
