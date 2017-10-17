/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ AFIndex, GNaIndex, mIndex, hIndex, jIndex, nVar};
static const char *compName = "INa";
static VARINFO varInfo[] =
{
   {"AF",PARAMETER_TYPE,AFIndex,0,0,1,1,"1"},
   {"GNa",PARAMETER_TYPE,GNaIndex,23,23,23,23,"mS/mF"},
   {"m",PSTATE_TYPE,mIndex,0,0,0,0,"1"},
   {"h",PSTATE_TYPE,hIndex,1,1,1,1,"1"},
   {"j",PSTATE_TYPE,jIndex,1,1,1,1,"1"}
};
typedef struct parameters_str { double  AF, GNa;} PARAMETERS;
typedef struct pstate_str { double  m, h, j;} PSTATE;
void Grandi_INaFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void Grandi_INaAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
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
         case GNaIndex:
            *value = parms->GNa; 
            break;
         case mIndex:
            *value = state->m; 
            break;
         case hIndex:
            *value = state->h; 
            break;
         case jIndex:
            *value = state->j; 
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
         case GNaIndex:
            parms->GNa = *value;
            break;
         case mIndex:
            state->m = *value;
            break;
         case hIndex:
            state->h = *value;
            break;
         case jIndex:
            state->j = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO Grandi_INaInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.pStateSize = sizeof(PSTATE);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 5;
   info.varInfo = varInfo;
   info.func = Grandi_INaFunc;
   info.access = Grandi_INaAccess;
   return info;
}
#ifdef doEnd
#define ENDCODE() {\
   if (derived->dState != 0) \
   {\
   double  *dState = derived->dState+pOffset;\
   dState[0]=dm; \
   dState[1]=dh; \
   dState[2]=dj; \
   }\
}
#else
#define ENDCODE()
#endif
