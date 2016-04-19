/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ GtoIndex, xtofIndex, ytofIndex, nVar};
static const char *compName = "Ito";
static VARINFO varInfo[] =
{
   {"Gto",PARAMETER_TYPE,GtoIndex,0.165,0.165,0.165,0.165,"mS/uF"},
   {"xtof",PSTATE_TYPE,xtofIndex,0.0,0.0,0.0,0.0,"1"},
   {"ytof",PSTATE_TYPE,ytofIndex,1.0,1.0,1.0,1.0,"1"}
};
typedef struct parameters_str { double  Gto;} PARAMETERS;
typedef struct pstate_str { double  xtof, ytof;} PSTATE;
void Grandi_ItoFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void Grandi_ItoAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PSTATE *state = (PSTATE *)statePtr;
   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case GtoIndex:
            *value = parms->Gto; 
            break;
         case xtofIndex:
            *value = state->xtof; 
            break;
         case ytofIndex:
            *value = state->ytof; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case GtoIndex:
            parms->Gto = *value;
            break;
         case xtofIndex:
            state->xtof = *value;
            break;
         case ytofIndex:
            state->ytof = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO Grandi_ItoInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.pStateSize = sizeof(PSTATE);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 3;
   info.varInfo = varInfo;
   info.func = Grandi_ItoFunc;
   info.access = Grandi_ItoAccess;
   return info;
}
#ifdef doEnd
#define ENDCODE() {\
   if (derived->dState != 0) \
   {\
   double  *dState = derived->dState+pOffset;\
   dState[0]=dxtof; \
   dState[1]=dytof; \
   }\
}
#else
#define ENDCODE()
#endif
