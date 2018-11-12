/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ GNaIndex, C3Index, C2Index, C1Index, OIndex, ISIndex, IC3Index, IC2Index, IFIndex, IM1Index, IM2Index, BC3Index, BC2Index, BC1Index, BOIndex, PaddingIndex, nVar};
static const char *compName = "INa";
static VARINFO varInfo[] =
{
   {"GNa",PARAMETER_TYPE,GNaIndex,13.5,13.5,13.5,"mS/uF"},
   {"C3",PSTATE_TYPE,C3Index,1,1,1,"1"},
   {"C2",PSTATE_TYPE,C2Index,0,0,0,"1"},
   {"C1",PSTATE_TYPE,C1Index,0,0,0,"1"},
   {"O",PSTATE_TYPE,OIndex,0,0,0,"1"},
   {"IS",PSTATE_TYPE,ISIndex,0,0,0,"1"},
   {"IC3",PSTATE_TYPE,IC3Index,0,0,0,"1"},
   {"IC2",PSTATE_TYPE,IC2Index,0,0,0,"1"},
   {"IF",PSTATE_TYPE,IFIndex,0,0,0,"1"},
   {"IM1",PSTATE_TYPE,IM1Index,0,0,0,"1"},
   {"IM2",PSTATE_TYPE,IM2Index,0,0,0,"1"},
   {"BC3",PSTATE_TYPE,BC3Index,0,0,0,"1"},
   {"BC2",PSTATE_TYPE,BC2Index,0,0,0,"1"},
   {"BC1",PSTATE_TYPE,BC1Index,0,0,0,"1"},
   {"BO",PSTATE_TYPE,BOIndex,0,0,0,"1"},
   {"Padding",PSTATE_TYPE,PaddingIndex,0,0,0,"1"}
};
typedef struct parameters_str { double  GNa;} PARAMETERS;
typedef struct pstate_str { double  C3, C2, C1, O, IS, IC3, IC2, IF, IM1, IM2, BC3, BC2, BC1, BO, Padding;} PSTATE;
void MYBGBKC_INaFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void MYBGBKC_INaAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PSTATE *state = (PSTATE *)statePtr;
   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case GNaIndex:
            *value = parms->GNa; 
            break;
         case C3Index:
            *value = state->C3; 
            break;
         case C2Index:
            *value = state->C2; 
            break;
         case C1Index:
            *value = state->C1; 
            break;
         case OIndex:
            *value = state->O; 
            break;
         case ISIndex:
            *value = state->IS; 
            break;
         case IC3Index:
            *value = state->IC3; 
            break;
         case IC2Index:
            *value = state->IC2; 
            break;
         case IFIndex:
            *value = state->IF; 
            break;
         case IM1Index:
            *value = state->IM1; 
            break;
         case IM2Index:
            *value = state->IM2; 
            break;
         case BC3Index:
            *value = state->BC3; 
            break;
         case BC2Index:
            *value = state->BC2; 
            break;
         case BC1Index:
            *value = state->BC1; 
            break;
         case BOIndex:
            *value = state->BO; 
            break;
         case PaddingIndex:
            *value = state->Padding; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case GNaIndex:
            parms->GNa = *value;
            break;
         case C3Index:
            state->C3 = *value;
            break;
         case C2Index:
            state->C2 = *value;
            break;
         case C1Index:
            state->C1 = *value;
            break;
         case OIndex:
            state->O = *value;
            break;
         case ISIndex:
            state->IS = *value;
            break;
         case IC3Index:
            state->IC3 = *value;
            break;
         case IC2Index:
            state->IC2 = *value;
            break;
         case IFIndex:
            state->IF = *value;
            break;
         case IM1Index:
            state->IM1 = *value;
            break;
         case IM2Index:
            state->IM2 = *value;
            break;
         case BC3Index:
            state->BC3 = *value;
            break;
         case BC2Index:
            state->BC2 = *value;
            break;
         case BC1Index:
            state->BC1 = *value;
            break;
         case BOIndex:
            state->BO = *value;
            break;
         case PaddingIndex:
            state->Padding = *value;
            break;
            assert(0); 
      }
   }
}
void   MYBGBKC_INaConstants();
COMPONENTINFO MYBGBKC_INaInit()
{
   COMPONENTINFO info;
   MYBGBKC_INaConstants();
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.pStateSize = sizeof(PSTATE);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 16;
   info.varInfo = varInfo;
   info.func = MYBGBKC_INaFunc;
   info.access = MYBGBKC_INaAccess;
   return info;
}
#ifdef doEnd
#define ENDCODE() {\
   if (derived->dState != 0) \
   {\
   double  *dState = derived->dState+pOffset;\
   dState[0]=dC3; \
   dState[1]=dC2; \
   dState[2]=dC1; \
   dState[3]=dO; \
   dState[4]=dIS; \
   dState[5]=dIC3; \
   dState[6]=dIC2; \
   dState[7]=dIF; \
   dState[8]=dIM1; \
   dState[9]=dIM2; \
   dState[10]=dBC3; \
   dState[11]=dBC2; \
   dState[12]=dBC1; \
   dState[13]=dBO; \
   dState[14]=dPadding; \
   }\
}
#else
#define ENDCODE()
#endif
