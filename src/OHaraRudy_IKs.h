/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ GKsIndex, Xs1Index, Xs2Index, nVar};
static const char *compName = "IKs";
static VARINFO varInfo[] =
{
   {"GKs",PARAMETER_TYPE,GKsIndex,0.0034,0.00476,0.0034,"mS/uF"},
   {"Xs1",PSTATE_TYPE,Xs1Index,0.270492,0.270492,0.270492,"1"},
   {"Xs2",PSTATE_TYPE,Xs2Index,0.0001963,0.0001963,0.0001963,"1"}
};
typedef struct parameters_str { double  GKs;} PARAMETERS;
typedef struct pstate_str { double  Xs1, Xs2;} PSTATE;
void OHaraRudy_IKsFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void OHaraRudy_IKsAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PSTATE *state = (PSTATE *)statePtr;
   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case GKsIndex:
            *value = parms->GKs; 
            break;
         case Xs1Index:
            *value = state->Xs1; 
            break;
         case Xs2Index:
            *value = state->Xs2; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case GKsIndex:
            parms->GKs = *value;
            break;
         case Xs1Index:
            state->Xs1 = *value;
            break;
         case Xs2Index:
            state->Xs2 = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO OHaraRudy_IKsInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.pStateSize = sizeof(PSTATE);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 3;
   info.varInfo = varInfo;
   info.func = OHaraRudy_IKsFunc;
   info.access = OHaraRudy_IKsAccess;
   return info;
}
#ifdef doEnd
#define ENDCODE() {\
   if (derived->dState != 0) \
   {\
   double  *dState = derived->dState+pOffset;\
   dState[0]=dXs1; \
   dState[1]=dXs2; \
   }\
}
#else
#define ENDCODE()
#endif
