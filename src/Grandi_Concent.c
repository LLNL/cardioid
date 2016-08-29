#include <math.h>
#include <stdlib.h>
#include "Grandi.h"
#include "Grandi_Concent.h"

void Grandi_ConcentFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt)
{
   VOLTAGE *voltage = (VOLTAGE *)state; 
   PSTATE *pState = (PSTATE *)(state+pOffset) ; 

   double Csqnb=pState->Csqnb; 
   double NaBj=pState->NaBj; 
   double NaBsl=pState->NaBsl; 
   double Naj=pState->Naj; 
   double Nasl=pState->Nasl; 
   double Nai=pState->Nai; 
   double Ki=pState->Ki; 
   double Casr=pState->Casr; 
   double Caj=pState->Caj; 
   double Casl=pState->Casl;
   double Cai=pState->Cai;

   //printf("Csqnb=%20.12f, NaBj=%20.12f, NaBsl=%20.12f, Naj=%20.12f, Nasl=%20.12f, Nai=%20.12f, Ki=%20.12f, Casr=%20.12f, Caj=%20.12f, Casl=%20.12f, Cai=%20.12f\n",Csqnb,NaBj,NaBsl,Naj,Nasl,Nai,Ki,Casr,Caj,Casl,Cai);
 
   CURRENTS *I = &(derived->I); 
   FLUXES J =derived->J; 
   I->stimulus = voltage->iStim; 

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
   double Bmax_Csqn=140.0e-3*Vmyo/Vsr;            // [mM] // Bmax_Csqn=2.6;      // Csqn buffering

   double INa_tot_junc=I->Na_junc+I->Nab_junc+3*I->NCX_junc+3*I->NaK_junc+I->CaNa_junc+I->NaL_junc;   // [uA/uF]
   double INa_tot_sl=I->Na_sl+I->Nab_sl+3*I->NCX_sl+3*I->NaK_sl+I->CaNa_sl+I->NaL_sl;   // [uA/uF]
   double ICa_tot_junc=I->Ca_junc+I->Cab_junc+I->pCa_junc-2*I->NCX_junc;                   // [uA/uF]
   double ICa_tot_sl=I->Ca_sl+I->Cab_sl+I->pCa_sl-2*I->NCX_sl;            // [uA/uF]

   double IK_tot=I->to+I->Kr+I->Ks_junc+I->Ks_sl+I->K1-2*(I->NaK_junc+I->NaK_sl)+I->CaK+I->Kp_junc+I->Kp_sl+I->Kur;     // [uA/uF] //SVP: added IKur
   double ICl_tot=I->ClCa_junc+I->ClCa_sl+I->Clb;                        // [uA/uF]
   double Itot=INa_tot_junc+INa_tot_sl+ICa_tot_junc+ICa_tot_sl+IK_tot+ICl_tot;
   double dVm=-Itot;

   //printf("ICab_junc=%20.12f\t ICab_sl=%20.12f\t\n",I->Cab_junc,I->Cab_sl);

   //printf("Itot=%20.12f\t INa_tot_junc=%20.12f\t INa_tot_sl=%20.12f\t ICa_tot_junc=%20.12f\t ICa_tot_sl=%20.12f\t IK_tot=%20.12f\t ICl_tot=%20.12f\t\n",Itot,INa_tot_junc,INa_tot_sl,ICa_tot_junc,ICa_tot_sl,IK_tot,ICl_tot);


   double dCsqnb=0;//kon_csqn*Casr*(Bmax_Csqn-Csqnb)-koff_csqn*Csqnb;       // Csqn      [mM/ms]
   double KmCsqnb=koff_csqn/kon_csqn;
   double Buff_Csqnb=1.0/(1.0+Bmax_Csqn*KmCsqnb/pow(KmCsqnb+Casr,2.0));

   //printf("Buff_Csqnb=%20.12f\n",Buff_Csqnb);

   double dNaBj=kon_na*Naj*(Bmax_Naj-NaBj)-koff_na*NaBj;        // NaBj      [mM/ms]
   double dNaBsl=kon_na*Nasl*(Bmax_Nasl-NaBsl)-koff_na*NaBsl;       // NaBsl     [mM/ms]
   double dNaj=-INa_tot_junc*Cmem/(Vjunc*F)+Jna_juncsl/Vjunc*(Nasl-Naj)-dNaBj;
   double dNasl=-INa_tot_sl*Cmem/(Vsl*F)+Jna_juncsl/Vsl*(Naj-Nasl)+Jna_slmyo/Vsl*(Nai-Nasl)-dNaBsl;
   double dNai=Jna_slmyo/Vmyo*(Nasl-Nai);             // [mM/msec]
   double dKi=-IK_tot*Cmem*(Vmyo*F);

   //double dCasr=J.up-(J.leak*Vmyo/Vsr+J.rel)-dCsqnb;         // Casr     [mM/ms] //Ratio 3 leak current
   double dCasr=Buff_Csqnb*( J.up-(J.leak*Vmyo/Vsr+J.rel) );

   double dCaj=-ICa_tot_junc*Cmem/(Vjunc*2*F)+Jca_juncsl/Vjunc*(Casl-Caj)-J.CaB_junc+(J.rel)*Vsr/Vjunc+J.leak*Vmyo/Vjunc;  // Ca_j
   double dCasl=-ICa_tot_sl*Cmem/(Vsl*2*F)+Jca_juncsl/Vsl*(Caj-Casl)+ Jca_slmyo/Vsl*(Cai-Casl)-J.CaB_sl;   // Ca_sl
   double dCai=-J.up*Vsr/Vmyo-J.CaB_cytosol +Jca_slmyo/Vmyo*(Casl-Cai);
   
   //printf("kon_csqn=%20.12f, Casr=%20.12f, Bmax_Csqn=%20.12f, Csqnb=%20.12f, koff_csqn=%20.12f\n",kon_csqn, Casr, Bmax_Csqn, Csqnb, koff_csqn);
   //printf("Vmyo/Vsr=%20.12f, Vsr=%20.12f, Vmyo=%20.12f\n",Vmyo/Vsr,Vsr,Vmyo);
   //printf("J.up=%20.12f, J.leak=%20.12f, J.rel=%20.12f, dCsqnb=%20.12f, ICa_tot_sl=%20.12f, Jca_juncsl=%20.12f, Jca_slmyo=%20.12f, J.CaB_sl\n\n",J.up,J.leak,J.rel,dCsqnb,ICa_tot_sl,Jca_juncsl,Jca_slmyo,J.CaB_sl);

   voltage->dVm = dVm; 
   derived->dState[0] = dVm; 

   ENDCODE()

     /*
   pState->Csqnb  *= exp(dt*dCsqnb/Csqnb);
   pState->NaBj  *= exp(dt*dNaBj/NaBj);
   pState->NaBsl  *= exp(dt*dNaBsl/NaBsl);
   pState->Naj  *= exp(dt*dNaj/Naj);
   pState->Nasl  *= exp(dt*dNasl/Nasl);
   pState->Nai  *= exp(dt*dNai/Nai);
   pState->Ki  *= exp(dt*dKi/Ki);
   pState->Casr  *= exp(dt*dCasr/Casr);
   pState->Caj  *= exp(dt*dCaj/Caj);
   pState->Casl  *= exp(dt*dCasl/Casl); 
   pState->Cai  *= exp(dt*dCai/Cai);
     */

   pState->Csqnb  += dt*dCsqnb;
   pState->NaBj  += dt*dNaBj;
   pState->NaBsl  += dt*dNaBsl;
   pState->Naj  += dt*dNaj;
   pState->Nasl  += dt*dNasl;
   pState->Nai  += dt*dNai;
   pState->Ki  += dt*dKi;
   pState->Casr  += dt*dCasr;
   pState->Caj  += dt*dCaj;
   pState->Casl  += dt*dCasl; 
   pState->Cai  += dt*dCai;
}
