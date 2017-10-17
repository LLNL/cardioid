/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ CaMKtrapIndex, nVar};
static const char *compName = "CaMKtrap";
static VARINFO varInfo[] =
{
   {"CaMKtrap",PSTATE_TYPE,CaMKtrapIndex,0.0124065,0.0124065,0.0124065,"1"}
};
typedef struct pstate_str { double  CaMKtrap;} PSTATE;
void OHaraRudy_CaMKtrapFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void OHaraRudy_CaMKtrapAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PSTATE *state = (PSTATE *)statePtr;
   if (type == READ)
   {
      switch (index)
      {
         case CaMKtrapIndex:
            *value = state->CaMKtrap; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case CaMKtrapIndex:
            state->CaMKtrap = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO OHaraRudy_CaMKtrapInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.pStateSize = sizeof(PSTATE);
   info.nVar = 1;
   info.varInfo = varInfo;
   info.func = OHaraRudy_CaMKtrapFunc;
   info.access = OHaraRudy_CaMKtrapAccess;
   return info;
}
#ifdef doEnd
#define ENDCODE() {\
   if (derived->dState != 0) \
   {\
   double  *dState = derived->dState+pOffset;\
   dState[0]=dCaMKtrap; \
   }\
}
#else
#define ENDCODE()
#endif
