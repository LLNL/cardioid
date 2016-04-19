#include <math.h>
#include "Grandi.h"
#include "Grandi_Fluxes.h"


void Grandi_FluxesFunc(CELLPARMS *parmsPtr, double *cell, int pOffset, DERIVED *derived, double dt )
{
   CONCENTRATIONS   *concentrations = (CONCENTRATIONS*) (cell + CONCENTRATIONS_OFFSET); 
   PSTATE *pState = (PSTATE *)(cell+pOffset) ; 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 

   double Naj = concentrations->Naj;
   double Nasl = concentrations->Nasl;
   double Casr = concentrations->Casr;
   double Caj = concentrations->Caj;
   double Casl = concentrations->Casl;
   double Cai = concentrations->Cai;
   double Ki = concentrations->Ki;

   double RyRr=pState->RyRr;
   double RyRo=pState->RyRo;
   double RyRi=pState->RyRi;
   double TnCL=pState->TnCL;
   double TnCHc=pState->TnCHc;
   double TnCHm=pState->TnCHm;
   double CaM=pState->CaM;
   double Myc=pState->Myc;
   double Mym=pState->Mym;
   double SRB=pState->SRB;
   double SLLj=pState->SLLj;
   double SLLsl=pState->SLLsl;
   double SLHj=pState->SLHj;
   double SLHsl=pState->SLHsl;

   //printf("RyRr=%20.12f, RyRo=%20.12f, RyRi=%20.12f, TnCL=%20.12f, TnCHc=%20.12f, TnCHm=%20.12f, CaM=%20.12f, Myc=%20.12f, Mym=%20.12f, SRB=%20.12f, SLLj=%20.12f, SLLsl=%20.12f, SLHj=%20.12f, SLHsl=%20.12f\n",RyRr,RyRo,RyRi,TnCL,TnCHc,TnCHm,CaM,Myc,Mym,SRB,SLLj,SLLsl,SLHj,SLHsl);

   //////Params not yet declared
   double Qpow=(T-310.0)/10.0;
   double Vcell=3.14159*cellRadius*cellRadius*cellLength*1.0e-15;    // [L]
   double Vmyo=0.65*Vcell;
   double Vsr=0.035*Vcell;
   double Vsl=0.02*Vcell;
   double Vjunc=0.0539*.01*Vcell;
   double Jca_juncsl =1.0/1.2134e12; // [L/msec]=8.2413e-13
   double Jca_slmyo=1.0/2.68510e11; // [L/msec]=3.2743e-12
   double Jna_juncsl=1.0/(1.6382e12/3.0*100.0); // [L/msec]=6.1043e-13
   double Jna_slmyo=1.0/(1.8308e10/3.0*100.0);  // [L/msec]=5.4621e-11

   ////// Jrel

   double Kmf=(2.5-1.25*ISO)*0.246e-3;          // [mM] default
   double koCa=10.0+20.0*AF+10.0*ISO*(1.0-AF);               // [mM^-2 1/ms]   //default 10   modified 20

   double kCaSR=MaxSR - (MaxSR-MinSR)/(1.0+pow(ec50SR/Casr,2.5));
   double koSRCa=koCa/kCaSR;
   double kiSRCa=kiCa*kCaSR;
   double RI=1-RyRr-RyRo-RyRi;

   derived->J.rel=cP->ks*RyRo*(Casr-Caj);          // [mM/ms]

   double dRyRr=(kim*RI-kiSRCa*Caj*RyRr)-(koSRCa*Caj*Caj*RyRr-kom*RyRo);   // R
   double dRyRo=(koSRCa*Caj*Caj*RyRr-kom*RyRo)-(kiSRCa*Caj*RyRo-kim*RyRi);// O
   double dRyRi=(kiSRCa*Caj*RyRo-kim*RyRi)-(kom*RyRi-koSRCa*Caj*Caj*RI);   // I

   ////// JSERCA

   derived->J.up=pow(Q10SRCaP,Qpow)*cP->Vmax_SRCaP*(pow(Cai/Kmf,hillSRCaP)-pow(Casr/Kmr,hillSRCaP))/(1.0+pow(Cai/Kmf,hillSRCaP)+pow(Casr/Kmr,hillSRCaP));

   ////// Jleak

   double phi=(1.0+0.25*AF);
   derived->J.leak=cP->Vmax_leak*phi*(Casr-Caj);           //   [mM/ms]

   ////// JBuffers

   double koff_tncl=(1.0+0.5*ISO)*19.6e-3;    // [1/ms]
   double Bmax_SR=19*.9e-3;     // [mM] (Bers text says 47e-3) 19e-3
   double Bmax_SLlowsl=37.4e-3*Vmyo/Vsl;        // [mM]    // SL buffering
   double Bmax_SLlowj=4.6e-3*Vmyo/Vjunc*0.1;    // [mM]    //Fei *0.1!!! junction reduction factor
   double Bmax_SLhighsl=13.4e-3*Vmyo/Vsl;       // [mM]
   double Bmax_SLhighj=1.65e-3*Vmyo/Vjunc*0.1;  // [mM] //Fei *0.1!!! junction reduction factor

   double dTnCL=kon_tncl*Cai*(Bmax_TnClow-TnCL)-koff_tncl*TnCL;            // TnCL      [mM/ms]
   double dTnCHc=kon_tnchca*Cai*(Bmax_TnChigh-TnCHc-TnCHm)-koff_tnchca*TnCHc; // TnCHc     [mM/ms]
   double dTnCHm=kon_tnchmg*Mgi*(Bmax_TnChigh-TnCHc-TnCHm)-koff_tnchmg*TnCHm;   // TnCHm     [mM/ms]
   double dCaM=kon_cam*Cai*(Bmax_CaM-CaM)-koff_cam*CaM;                 // CaM       [mM/ms]
   double dMyc=kon_myoca*Cai*(Bmax_myosin-Myc-Mym)-koff_myoca*Myc;    // Myosin_ca [mM/ms]
   double dMym=kon_myomg*Mgi*(Bmax_myosin-Myc-Mym)-koff_myomg*Mym;      // Myosin_mg [mM/ms]
   double dSRB=kon_sr*Cai*(Bmax_SR-SRB)-koff_sr*SRB;                    // SRB       [mM/ms]
   double dSLLj=kon_sll*Caj*(Bmax_SLlowj-SLLj)-koff_sll*SLLj;       // SLLj      [mM/ms]
   double dSLLsl=kon_sll*Casl*(Bmax_SLlowsl-SLLsl)-koff_sll*SLLsl;      // SLLsl     [mM/ms]
   double dSLHj=kon_slh*Caj*(Bmax_SLhighj-SLHj)-koff_slh*SLHj;      // SLHj      [mM/ms]
   double dSLHsl=kon_slh*Casl*(Bmax_SLhighsl-SLHsl)-koff_slh*SLHsl;     // SLHsl     [mM/ms]

   derived->J.CaB_cytosol=dTnCL+dTnCHc+dTnCHm+dCaM+dMyc+dMym+dSRB;
   derived->J.CaB_junc=dSLLj+dSLHj;
   derived->J.CaB_sl=dSLLsl+dSLHsl;

   //   printf("J.up=%20.12f, J.leak=%20.12f, J.rel=%20.12f, J.CaB_junc=%20.12f, J.CaB_sl=%20.12f, J.CaB_cytosol=%20.12f\n",derived->J.up,derived->J.leak,derived->J.rel,derived->J.CaB_junc,derived->J.CaB_sl,derived->J.CaB_cytosol);

   //////
   
   ENDCODE()

   pState->RyRr += dt*dRyRr;
   pState->RyRo += dt*dRyRo;
   pState->RyRi += dt*dRyRi;
   pState->TnCL += dt*dTnCL;
   pState->TnCHc += dt*dTnCHc;
   pState->TnCHm += dt*dTnCHm;
   pState->CaM += dt*dCaM;
   pState->Myc += dt*dMyc;
   pState->Mym += dt*dMym;
   pState->SRB += dt*dSRB;
   pState->SLLj += dt*dSLLj;
   pState->SLLsl += dt*dSLLsl;
   pState->SLHj += dt*dSLHj;
   pState->SLHsl += dt*dSLHsl;

/*      pState->RyRr*=exp(dt*dRyRr/RyRr); */
/*      pState->RyRo*=exp(dt*dRyRo/RyRo); */
/*      pState->RyRi*=exp(dt*dRyRi/RyRi); */
/*      pState->TnCL*=exp(dt*dTnCL/TnCL); */
/*      pState->TnCHc*=exp(dt*dTnCHc/TnCHc); */
/*      pState->TnCHm*=exp(dt*dTnCHm/TnCHm); */
/*      pState->CaM*=exp(dt*dCaM/CaM); */
/*      pState->Myc*=exp(dt*dMyc/Myc); */
/*      pState->Mym*=exp(dt*dMym/Mym); */
/*      pState->SRB*=exp(dt*dSRB/SRB); */
/*      pState->SLLj*=exp(dt*dSLLj/SLLj); */
/*      pState->SLLsl*=exp(dt*dSLLsl/SLLsl); */
/*      pState->SLHj*=exp(dt*dSLHj/SLHj); */
/*      pState->SLHsl*=exp(dt*dSLHsl/SLHsl); */
}
