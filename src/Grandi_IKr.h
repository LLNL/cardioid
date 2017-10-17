/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ GKrIndex, xkrIndex, nVar};
static const char *compName = "IKr";
static VARINFO varInfo[] =
{
   {"GKr",PARAMETER_TYPE,GKrIndex,0.035,0.035,0.035,0.035,"mS/uF"},
   {"xkr",PSTATE_TYPE,xkrIndex,0,0,0,0,"1"}
};
typedef struct parameters_str { double  GKr;} PARAMETERS;
typedef struct pstate_str { double  xkr;} PSTATE;
void Grandi_IKrFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void Grandi_IKrAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
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
         case xkrIndex:
            *value = state->xkr; 
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
         case xkrIndex:
            state->xkr = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO Grandi_IKrInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.pStateSize = sizeof(PSTATE);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 2;
   info.varInfo = varInfo;
   info.func = Grandi_IKrFunc;
   info.access = Grandi_IKrAccess;
   return info;
}
#ifdef doEnd
#define ENDCODE() {\
   if (derived->dState != 0) \
   {\
   double  *dState = derived->dState+pOffset;\
   dState[0]=dxkr; \
   }\
}
#else
#define ENDCODE()
#endif
