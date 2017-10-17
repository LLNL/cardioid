/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ AFIndex, GKsIndex, xksIndex, nVar};
static const char *compName = "IKs";
static VARINFO varInfo[] =
{
   {"AF",PARAMETER_TYPE,AFIndex,0,0,1,1,"1"},
   {"GKs",PARAMETER_TYPE,GKsIndex,0.0035,0.0035,0.0035,0.0035,"mS/uF"},
   {"xks",PSTATE_TYPE,xksIndex,0,0,0,0,"1"}
};
typedef struct parameters_str { double  AF, GKs;} PARAMETERS;
typedef struct pstate_str { double  xks;} PSTATE;
void Grandi_IKsFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void Grandi_IKsAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
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
         case GKsIndex:
            *value = parms->GKs; 
            break;
         case xksIndex:
            *value = state->xks; 
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
         case GKsIndex:
            parms->GKs = *value;
            break;
         case xksIndex:
            state->xks = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO Grandi_IKsInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.pStateSize = sizeof(PSTATE);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 3;
   info.varInfo = varInfo;
   info.func = Grandi_IKsFunc;
   info.access = Grandi_IKsAccess;
   return info;
}
#ifdef doEnd
#define ENDCODE() {\
   if (derived->dState != 0) \
   {\
   double  *dState = derived->dState+pOffset;\
   dState[0]=dxks; \
   }\
}
#else
#define ENDCODE()
#endif
