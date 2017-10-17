/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ CsqnbIndex, NaBjIndex, NaBslIndex, NajIndex, NaslIndex, NaiIndex, KiIndex, CasrIndex, CajIndex, CaslIndex, CaiIndex, nVar};
static const char *compName = "Concent";
static VARINFO varInfo[] =
{
   {"Csqnb",PSTATE_TYPE,CsqnbIndex,1.25,1.25,1.25,1.25,"mM"},
   {"NaBj",PSTATE_TYPE,NaBjIndex,3.5,3.5,3.5,3.5,"mM"},
   {"NaBsl",PSTATE_TYPE,NaBslIndex,0.8,0.8,0.8,0.8,"mM"},
   {"Naj",PSTATE_TYPE,NajIndex,9.136,9.136,9.136,9.136,"mM"},
   {"Nasl",PSTATE_TYPE,NaslIndex,9.136,9.136,9.136,9.136,"mM"},
   {"Nai",PSTATE_TYPE,NaiIndex,9.136,9.136,9.136,9.136,"mM"},
   {"Ki",PSTATE_TYPE,KiIndex,120,120,120,120,"mM"},
   {"Casr",PSTATE_TYPE,CasrIndex,0.01,0.01,0.01,0.01,"mM"},
   {"Caj",PSTATE_TYPE,CajIndex,0.00017,0.00017,0.00017,0.00017,"mM"},
   {"Casl",PSTATE_TYPE,CaslIndex,0.0001,0.0001,0.0001,0.0001,"mM"},
   {"Cai",PSTATE_TYPE,CaiIndex,0.0001,0.0001,0.0001,0.0001,"mM"}
};
typedef struct pstate_str { double  Csqnb, NaBj, NaBsl, Naj, Nasl, Nai, Ki, Casr, Caj, Casl, Cai;} PSTATE;
void Grandi_ConcentFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void Grandi_ConcentAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PSTATE *state = (PSTATE *)statePtr;
   if (type == READ)
   {
      switch (index)
      {
         case CsqnbIndex:
            *value = state->Csqnb; 
            break;
         case NaBjIndex:
            *value = state->NaBj; 
            break;
         case NaBslIndex:
            *value = state->NaBsl; 
            break;
         case NajIndex:
            *value = state->Naj; 
            break;
         case NaslIndex:
            *value = state->Nasl; 
            break;
         case NaiIndex:
            *value = state->Nai; 
            break;
         case KiIndex:
            *value = state->Ki; 
            break;
         case CasrIndex:
            *value = state->Casr; 
            break;
         case CajIndex:
            *value = state->Caj; 
            break;
         case CaslIndex:
            *value = state->Casl; 
            break;
         case CaiIndex:
            *value = state->Cai; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case CsqnbIndex:
            state->Csqnb = *value;
            break;
         case NaBjIndex:
            state->NaBj = *value;
            break;
         case NaBslIndex:
            state->NaBsl = *value;
            break;
         case NajIndex:
            state->Naj = *value;
            break;
         case NaslIndex:
            state->Nasl = *value;
            break;
         case NaiIndex:
            state->Nai = *value;
            break;
         case KiIndex:
            state->Ki = *value;
            break;
         case CasrIndex:
            state->Casr = *value;
            break;
         case CajIndex:
            state->Caj = *value;
            break;
         case CaslIndex:
            state->Casl = *value;
            break;
         case CaiIndex:
            state->Cai = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO Grandi_ConcentInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.pStateSize = sizeof(PSTATE);
   info.nVar = 11;
   info.varInfo = varInfo;
   info.func = Grandi_ConcentFunc;
   info.access = Grandi_ConcentAccess;
   return info;
}
#ifdef doEnd
#define ENDCODE() {\
   if (derived->dState != 0) \
   {\
   double  *dState = derived->dState+pOffset;\
   dState[0]=dCsqnb; \
   dState[1]=dNaBj; \
   dState[2]=dNaBsl; \
   dState[3]=dNaj; \
   dState[4]=dNasl; \
   dState[5]=dNai; \
   dState[6]=dKi; \
   dState[7]=dCasr; \
   dState[8]=dCaj; \
   dState[9]=dCasl; \
   dState[10]=dCai; \
   }\
}
#else
#define ENDCODE()
#endif
