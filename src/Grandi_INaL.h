/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ AFIndex, GNaLIndex, mLIndex, hLIndex, nVar};
static const char *compName = "INaL";
static VARINFO varInfo[] =
{
   {"AF",PARAMETER_TYPE,AFIndex,0,0,1,1,"1"},
   {"GNaL",PARAMETER_TYPE,GNaLIndex,0.0025,0.0025,0.0025,0.0025,"mS/mF"},
   {"mL",PSTATE_TYPE,mLIndex,0,0,0,0,"1"},
   {"hL",PSTATE_TYPE,hLIndex,1,1,1,1,"1"}
};
typedef struct parameters_str { double  AF, GNaL;} PARAMETERS;
typedef struct pstate_str { double  mL, hL;} PSTATE;
void Grandi_INaLFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void Grandi_INaLAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
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
         case GNaLIndex:
            *value = parms->GNaL; 
            break;
         case mLIndex:
            *value = state->mL; 
            break;
         case hLIndex:
            *value = state->hL; 
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
         case GNaLIndex:
            parms->GNaL = *value;
            break;
         case mLIndex:
            state->mL = *value;
            break;
         case hLIndex:
            state->hL = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO Grandi_INaLInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.pStateSize = sizeof(PSTATE);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 4;
   info.varInfo = varInfo;
   info.func = Grandi_INaLFunc;
   info.access = Grandi_INaLAccess;
   return info;
}
#ifdef doEnd
#define ENDCODE() {\
   if (derived->dState != 0) \
   {\
   double  *dState = derived->dState+pOffset;\
   dState[0]=dmL; \
   dState[1]=dhL; \
   }\
}
#else
#define ENDCODE()
#endif
