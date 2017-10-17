/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ GK1Index, XK1Index, nVar};
static const char *compName = "IK1";
static VARINFO varInfo[] =
{
   {"GK1",PARAMETER_TYPE,GK1Index,0.1908,0.22896,0.24804,"mS/uF"},
   {"XK1",PSTATE_TYPE,XK1Index,0.996801,0.996801,0.996801,"1"}
};
typedef struct parameters_str { double  GK1;} PARAMETERS;
typedef struct pstate_str { double  XK1;} PSTATE;
void OHaraRudy_IK1Func(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void OHaraRudy_IK1Access(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PSTATE *state = (PSTATE *)statePtr;
   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case GK1Index:
            *value = parms->GK1; 
            break;
         case XK1Index:
            *value = state->XK1; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case GK1Index:
            parms->GK1 = *value;
            break;
         case XK1Index:
            state->XK1 = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO OHaraRudy_IK1Init()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.pStateSize = sizeof(PSTATE);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 2;
   info.varInfo = varInfo;
   info.func = OHaraRudy_IK1Func;
   info.access = OHaraRudy_IK1Access;
   return info;
}
#ifdef doEnd
#define ENDCODE() {\
   if (derived->dState != 0) \
   {\
   double  *dState = derived->dState+pOffset;\
   dState[0]=dXK1; \
   }\
}
#else
#define ENDCODE()
#endif
