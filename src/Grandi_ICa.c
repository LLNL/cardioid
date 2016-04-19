#include <math.h>
#include <assert.h>
#include "Grandi.h"
#include "Grandi_ICa.h"
void Grandi_ICaFunc(CELLPARMS *parmsPtr, double  *cell, int pOffset, DERIVED *derived, double dt )
{
   VOLTAGE *voltage = (VOLTAGE *)cell; 
   CONCENTRATIONS  *concentrations = (CONCENTRATIONS*) (cell + CONCENTRATIONS_OFFSET); 
   double v = voltage->Vm; 
   double Caj=concentrations->Caj;
   double Casl=concentrations->Casl;
   double Naj=concentrations->Naj;
   double Nasl=concentrations->Nasl;
   double Ki=concentrations->Ki;

   double phi=(1.0+0.5*ISO)*(1.0-0.5*AF);
   double Fsl_CaL=1.0-Fjunc_CaL;

   PSTATE *pState = (PSTATE *)(cell+pOffset) ; 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   
   double d=pState->d;
   double f=pState->f;
   double fcaBj=pState->fcaBj;
   double fcaBsl=pState->fcaBsl;

   double VFRT = v*FRT; 
   double VFFRT = F*VFRT;
   double exVFRT = exp(VFRT);
   double ex2VFRT = exp(2.0*VFRT);
   
   double ibarca_j=cP->pCa*4.0*VFFRT*(0.341*Caj*ex2VFRT-0.341*Cao)/(ex2VFRT-1.0);
   double ibarca_sl=cP->pCa*4.0*VFFRT*(0.341*Casl*ex2VFRT-0.341*Cao)/(ex2VFRT-1.0);
   double ibark=cP->pK*VFFRT*(0.75*Ki*exVFRT-0.75*Ko) /(exVFRT-1.0);
   double ibarna_j=cP->pNa*VFFRT*(0.75*Naj*exVFRT-0.75*Nao)/(exVFRT-1.0);
   double ibarna_sl=cP->pNa*VFFRT*(0.75*Nasl*exVFRT-0.75*Nao)/(exVFRT-1.0);

   double Q10CaL=1.8;
   double Qpow=(T-310.0)/10.0;

   derived->I.Ca_junc=(Fjunc_CaL*ibarca_j*d*f*(1-fcaBj)*pow(Q10CaL,Qpow))*0.45;
   derived->I.Ca_sl=(Fsl_CaL*ibarca_sl*d*f*(1-fcaBsl)*pow(Q10CaL,Qpow))*0.45;
   derived->I.CaNa_junc=(Fjunc_CaL*ibarna_j*d*f*(1-fcaBj)*pow(Q10CaL,Qpow))*0.45;
   derived->I.CaNa_sl=(Fsl_CaL*ibarna_sl*d*f*(1-fcaBsl)*pow(Q10CaL,Qpow))*0.45;
   derived->I.CaK=(ibark*d*f*(Fjunc_CaL*(1-fcaBj)+Fsl_CaL*(1-fcaBsl))*pow(Q10CaL,Qpow))*0.45;

   double dss=1.0/(1.0+exp(-(v+3.0*ISO+9.0)/6.0));
   double taud=1.0*dss*(1.0-exp(-(v+3.0*ISO+9.0)/6.0))/(0.035*(v+3.0*ISO+9.0));
   double fss=1.0/(1.0+exp((v+3.0*ISO+30.0)/7.0))+0.2/(1+exp((50.0-v-3.0*ISO)/20.0));
   double tauf=1.0/(0.0197*exp(-pow(0.0337*(v+3.0*ISO+25.0),2.0))+0.02);

   double dd=(dss-d)/taud;
   double df=(fss-f)/tauf;
   double dfcaBj=1.7*Caj*(1.0-fcaBj)-11.9e-3*fcaBj;
   double dfcaBsl=1.7*Casl*(1-fcaBsl)-11.9e-3*fcaBsl;

   ENDCODE()
     //pState->d += dt*dd;
     //pState->f += dt*df;
   pState->d=dss-(dss-d)*exp(-dt/taud);
   pState->f=fss-(fss-f)*exp(-dt/tauf);

   pState->fcaBj += dt*dfcaBj;
   pState->fcaBsl += dt*dfcaBsl;

}
