/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ GtoIndex, aDeltaIndex, aIndex, iFastIndex, iSlowIndex, aCaMKIndex, iCaMKFastIndex, iCaMKSlowIndex, nVar};
static const char *compName = "Ito";
static VARINFO varInfo[] =
{
   {"Gto",PARAMETER_TYPE,GtoIndex,0.02,0.08,0.08,"mS/uF"},
   {"aDelta",PARAMETER_TYPE,aDeltaIndex,0,0.95,0,"1"},
   {"a",PSTATE_TYPE,aIndex,0.00101185,0.00101185,0.00101185,"1"},
   {"iFast",PSTATE_TYPE,iFastIndex,0.999542,0.999542,0.999542,"1"},
   {"iSlow",PSTATE_TYPE,iSlowIndex,0.589579,0.589579,0.589579,"1"},
   {"aCaMK",PSTATE_TYPE,aCaMKIndex,0.000515567,0.000515567,0.000515567,"1"},
   {"iCaMKFast",PSTATE_TYPE,iCaMKFastIndex,0.999542,0.999542,0.999542,"1"},
   {"iCaMKSlow",PSTATE_TYPE,iCaMKSlowIndex,0.641861,0.641861,0.641861,"1"}
};
typedef struct parameters_str { double  Gto, aDelta;} PARAMETERS;
typedef struct pstate_str { double  a, iFast, iSlow, aCaMK, iCaMKFast, iCaMKSlow;} PSTATE;
void OHaraRudy_ItoFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void OHaraRudy_ItoAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
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
         case aDeltaIndex:
            *value = parms->aDelta; 
            break;
         case aIndex:
            *value = state->a; 
            break;
         case iFastIndex:
            *value = state->iFast; 
            break;
         case iSlowIndex:
            *value = state->iSlow; 
            break;
         case aCaMKIndex:
            *value = state->aCaMK; 
            break;
         case iCaMKFastIndex:
            *value = state->iCaMKFast; 
            break;
         case iCaMKSlowIndex:
            *value = state->iCaMKSlow; 
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
         case aDeltaIndex:
            parms->aDelta = *value;
            break;
         case aIndex:
            state->a = *value;
            break;
         case iFastIndex:
            state->iFast = *value;
            break;
         case iSlowIndex:
            state->iSlow = *value;
            break;
         case aCaMKIndex:
            state->aCaMK = *value;
            break;
         case iCaMKFastIndex:
            state->iCaMKFast = *value;
            break;
         case iCaMKSlowIndex:
            state->iCaMKSlow = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO OHaraRudy_ItoInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.pStateSize = sizeof(PSTATE);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 8;
   info.varInfo = varInfo;
   info.func = OHaraRudy_ItoFunc;
   info.access = OHaraRudy_ItoAccess;
   return info;
}
#ifdef doEnd
#define ENDCODE() {\
   if (derived->dState != 0) \
   {\
   double  *dState = derived->dState+pOffset;\
   dState[0]=da; \
   dState[1]=diFast; \
   dState[2]=diSlow; \
   dState[3]=daCaMK; \
   dState[4]=diCaMKFast; \
   dState[5]=diCaMKSlow; \
   }\
}
#else
#define ENDCODE()
#endif
