/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ AFIndex, ksIndex, Vmax_SRCaPIndex, Vmax_leakIndex, RyRrIndex, RyRoIndex, RyRiIndex, TnCLIndex, TnCHcIndex, TnCHmIndex, CaMIndex, MycIndex, MymIndex, SRBIndex, SLLjIndex, SLLslIndex, SLHjIndex, SLHslIndex, nVar};
static const char *compName = "Fluxes";
static VARINFO varInfo[] =
{
   {"AF",PARAMETER_TYPE,AFIndex,0,0,1,1,"1"},
   {"ks",PARAMETER_TYPE,ksIndex,25,25,25,25,"1/ms"},
   {"Vmax_SRCaP",PARAMETER_TYPE,Vmax_SRCaPIndex,0.0053114,0.0053114,0.0053114,0.0053114,"mM/ms"},
   {"Vmax_leak",PARAMETER_TYPE,Vmax_leakIndex,5.348e-06,5.348e-06,5.348e-06,5.348e-06,"mM/ms"},
   {"RyRr",PSTATE_TYPE,RyRrIndex,1,1,1,1,"1"},
   {"RyRo",PSTATE_TYPE,RyRoIndex,0,0,0,0,"1"},
   {"RyRi",PSTATE_TYPE,RyRiIndex,0,0,0,0,"1"},
   {"TnCL",PSTATE_TYPE,TnCLIndex,0.01,0.01,0.01,0.01,"1"},
   {"TnCHc",PSTATE_TYPE,TnCHcIndex,0.1,0.1,0.1,0.1,"1"},
   {"TnCHm",PSTATE_TYPE,TnCHmIndex,0.01,0.01,0.01,0.01,"1"},
   {"CaM",PSTATE_TYPE,CaMIndex,0.0003,0.0003,0.0003,0.0003,"1"},
   {"Myc",PSTATE_TYPE,MycIndex,0.0013,0.0013,0.0013,0.0013,"1"},
   {"Mym",PSTATE_TYPE,MymIndex,0.14,0.14,0.14,0.14,"1"},
   {"SRB",PSTATE_TYPE,SRBIndex,0.002,0.002,0.002,0.002,"1"},
   {"SLLj",PSTATE_TYPE,SLLjIndex,0.01,0.01,0.01,0.01,"1"},
   {"SLLsl",PSTATE_TYPE,SLLslIndex,0.1,0.1,0.1,0.1,"1"},
   {"SLHj",PSTATE_TYPE,SLHjIndex,0.0073,0.0073,0.0073,0.0073,"1"},
   {"SLHsl",PSTATE_TYPE,SLHslIndex,0.073,0.073,0.073,0.073,"1"}
};
typedef struct parameters_str { double  AF, ks, Vmax_SRCaP, Vmax_leak;} PARAMETERS;
typedef struct pstate_str { double  RyRr, RyRo, RyRi, TnCL, TnCHc, TnCHm, CaM, Myc, Mym, SRB, SLLj, SLLsl, SLHj, SLHsl;} PSTATE;
void Grandi_FluxesFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void Grandi_FluxesAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PSTATE *state = (PSTATE *)statePtr;
   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case AFIndex:
            *value = parms->AF; 
            break;
         case ksIndex:
            *value = parms->ks; 
            break;
         case Vmax_SRCaPIndex:
            *value = parms->Vmax_SRCaP; 
            break;
         case Vmax_leakIndex:
            *value = parms->Vmax_leak; 
            break;
         case RyRrIndex:
            *value = state->RyRr; 
            break;
         case RyRoIndex:
            *value = state->RyRo; 
            break;
         case RyRiIndex:
            *value = state->RyRi; 
            break;
         case TnCLIndex:
            *value = state->TnCL; 
            break;
         case TnCHcIndex:
            *value = state->TnCHc; 
            break;
         case TnCHmIndex:
            *value = state->TnCHm; 
            break;
         case CaMIndex:
            *value = state->CaM; 
            break;
         case MycIndex:
            *value = state->Myc; 
            break;
         case MymIndex:
            *value = state->Mym; 
            break;
         case SRBIndex:
            *value = state->SRB; 
            break;
         case SLLjIndex:
            *value = state->SLLj; 
            break;
         case SLLslIndex:
            *value = state->SLLsl; 
            break;
         case SLHjIndex:
            *value = state->SLHj; 
            break;
         case SLHslIndex:
            *value = state->SLHsl; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case AFIndex:
            parms->AF = *value;
            break;
         case ksIndex:
            parms->ks = *value;
            break;
         case Vmax_SRCaPIndex:
            parms->Vmax_SRCaP = *value;
            break;
         case Vmax_leakIndex:
            parms->Vmax_leak = *value;
            break;
         case RyRrIndex:
            state->RyRr = *value;
            break;
         case RyRoIndex:
            state->RyRo = *value;
            break;
         case RyRiIndex:
            state->RyRi = *value;
            break;
         case TnCLIndex:
            state->TnCL = *value;
            break;
         case TnCHcIndex:
            state->TnCHc = *value;
            break;
         case TnCHmIndex:
            state->TnCHm = *value;
            break;
         case CaMIndex:
            state->CaM = *value;
            break;
         case MycIndex:
            state->Myc = *value;
            break;
         case MymIndex:
            state->Mym = *value;
            break;
         case SRBIndex:
            state->SRB = *value;
            break;
         case SLLjIndex:
            state->SLLj = *value;
            break;
         case SLLslIndex:
            state->SLLsl = *value;
            break;
         case SLHjIndex:
            state->SLHj = *value;
            break;
         case SLHslIndex:
            state->SLHsl = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO Grandi_FluxesInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.pStateSize = sizeof(PSTATE);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 18;
   info.varInfo = varInfo;
   info.func = Grandi_FluxesFunc;
   info.access = Grandi_FluxesAccess;
   return info;
}
#ifdef doEnd
#define ENDCODE() {\
   if (derived->dState != 0) \
   {\
   double  *dState = derived->dState+pOffset;\
   dState[0]=dRyRr; \
   dState[1]=dRyRo; \
   dState[2]=dRyRi; \
   dState[3]=dTnCL; \
   dState[4]=dTnCHc; \
   dState[5]=dTnCHm; \
   dState[6]=dCaM; \
   dState[7]=dMyc; \
   dState[8]=dMym; \
   dState[9]=dSRB; \
   dState[10]=dSLLj; \
   dState[11]=dSLLsl; \
   dState[12]=dSLHj; \
   dState[13]=dSLHsl; \
   }\
}
#else
#define ENDCODE()
#endif
