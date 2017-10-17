/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ GKrIndex, DIndex, kCIndex, kOIndex, kIIndex, rCIndex, rOIndex, rIIndex, C3Index, C2Index, C1Index, OIndex, IIndex, C3dIndex, C2dIndex, C1dIndex, OdIndex, IdIndex, nVar};
static const char *compName = "IKr";
static VARINFO varInfo[] =
{
   {"GKr",PARAMETER_TYPE,GKrIndex,0.0552,0.07176,0.04416,"mS/uF"},
   {"D",PARAMETER_TYPE,DIndex,0,0,0,"mM"},
   {"kC",PARAMETER_TYPE,kCIndex,0,0,0,"1/(mM*ms)"},
   {"kO",PARAMETER_TYPE,kOIndex,0,0,0,"1/(mM*ms)"},
   {"kI",PARAMETER_TYPE,kIIndex,0,0,0,"1/(mM*ms)"},
   {"rC",PARAMETER_TYPE,rCIndex,0,0,0,"1/ms"},
   {"rO",PARAMETER_TYPE,rOIndex,0,0,0,"1/ms"},
   {"rI",PARAMETER_TYPE,rIIndex,0,0,0,"1/ms"},
   {"C3",PSTATE_TYPE,C3Index,1,1,1,"1"},
   {"C2",PSTATE_TYPE,C2Index,0,0,0,"1"},
   {"C1",PSTATE_TYPE,C1Index,0,0,0,"1"},
   {"O",PSTATE_TYPE,OIndex,0,0,0,"1"},
   {"I",PSTATE_TYPE,IIndex,0,0,0,"1"},
   {"C3d",PSTATE_TYPE,C3dIndex,0,0,0,"1"},
   {"C2d",PSTATE_TYPE,C2dIndex,0,0,0,"1"},
   {"C1d",PSTATE_TYPE,C1dIndex,0,0,0,"1"},
   {"Od",PSTATE_TYPE,OdIndex,0,0,0,"1"},
   {"Id",PSTATE_TYPE,IdIndex,0,0,0,"1"}
};
typedef struct parameters_str { double  GKr, D, kC, kO, kI, rC, rO, rI;} PARAMETERS;
typedef struct pstate_str { double  C3, C2, C1, O, I, C3d, C2d, C1d, Od, Id;} PSTATE;
void RTYSC14A_IKrFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void RTYSC14A_IKrAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PSTATE *state = (PSTATE *)statePtr;
   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case GKrIndex:
            *value = parms->GKr; 
            break;
         case DIndex:
            *value = parms->D; 
            break;
         case kCIndex:
            *value = parms->kC; 
            break;
         case kOIndex:
            *value = parms->kO; 
            break;
         case kIIndex:
            *value = parms->kI; 
            break;
         case rCIndex:
            *value = parms->rC; 
            break;
         case rOIndex:
            *value = parms->rO; 
            break;
         case rIIndex:
            *value = parms->rI; 
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
         case IIndex:
            *value = state->I; 
            break;
         case C3dIndex:
            *value = state->C3d; 
            break;
         case C2dIndex:
            *value = state->C2d; 
            break;
         case C1dIndex:
            *value = state->C1d; 
            break;
         case OdIndex:
            *value = state->Od; 
            break;
         case IdIndex:
            *value = state->Id; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case GKrIndex:
            parms->GKr = *value;
            break;
         case DIndex:
            parms->D = *value;
            break;
         case kCIndex:
            parms->kC = *value;
            break;
         case kOIndex:
            parms->kO = *value;
            break;
         case kIIndex:
            parms->kI = *value;
            break;
         case rCIndex:
            parms->rC = *value;
            break;
         case rOIndex:
            parms->rO = *value;
            break;
         case rIIndex:
            parms->rI = *value;
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
         case IIndex:
            state->I = *value;
            break;
         case C3dIndex:
            state->C3d = *value;
            break;
         case C2dIndex:
            state->C2d = *value;
            break;
         case C1dIndex:
            state->C1d = *value;
            break;
         case OdIndex:
            state->Od = *value;
            break;
         case IdIndex:
            state->Id = *value;
            break;
            assert(0); 
      }
   }
}
void   RTYSC14A_IKrConstants();
COMPONENTINFO RTYSC14A_IKrInit()
{
   COMPONENTINFO info;
   RTYSC14A_IKrConstants();
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.pStateSize = sizeof(PSTATE);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 18;
   info.varInfo = varInfo;
   info.func = RTYSC14A_IKrFunc;
   info.access = RTYSC14A_IKrAccess;
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
   dState[4]=dI; \
   dState[5]=dC3d; \
   dState[6]=dC2d; \
   dState[7]=dC1d; \
   dState[8]=dOd; \
   dState[9]=dId; \
   }\
}
#else
#define ENDCODE()
#endif
