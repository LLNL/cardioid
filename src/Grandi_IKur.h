/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ RAIndex, AFIndex, GKurIndex, xkurIndex, ykurIndex, nVar};
static const char *compName = "IKur";
static VARINFO varInfo[] =
{
   {"RA",PARAMETER_TYPE,RAIndex,1,0,1,0,"1"},
   {"AF",PARAMETER_TYPE,AFIndex,0,0,1,1,"1"},
   {"GKur",PARAMETER_TYPE,GKurIndex,0.045,0.045,0.045,0.045,"mS/uF"},
   {"xkur",PSTATE_TYPE,xkurIndex,0,0,0,0,"1"},
   {"ykur",PSTATE_TYPE,ykurIndex,1,1,1,1,"1"}
};
typedef struct parameters_str { double  RA, AF, GKur;} PARAMETERS;
typedef struct pstate_str { double  xkur, ykur;} PSTATE;
void Grandi_IKurFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void Grandi_IKurAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PSTATE *state = (PSTATE *)statePtr;
   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case RAIndex:
            *value = parms->RA; 
            break;
         case AFIndex:
            *value = parms->AF; 
            break;
         case GKurIndex:
            *value = parms->GKur; 
            break;
         case xkurIndex:
            *value = state->xkur; 
            break;
         case ykurIndex:
            *value = state->ykur; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case RAIndex:
            parms->RA = *value;
            break;
         case AFIndex:
            parms->AF = *value;
            break;
         case GKurIndex:
            parms->GKur = *value;
            break;
         case xkurIndex:
            state->xkur = *value;
            break;
         case ykurIndex:
            state->ykur = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO Grandi_IKurInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.pStateSize = sizeof(PSTATE);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 5;
   info.varInfo = varInfo;
   info.func = Grandi_IKurFunc;
   info.access = Grandi_IKurAccess;
   return info;
}
#ifdef doEnd
#define ENDCODE() {\
   if (derived->dState != 0) \
   {\
   double  *dState = derived->dState+pOffset;\
   dState[0]=dxkur; \
   dState[1]=dykur; \
   }\
}
#else
#define ENDCODE()
#endif
