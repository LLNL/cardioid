/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ GKbIndex, nVar};
static const char *compName = "IKb";
static VARINFO varInfo[] =
{
   {"GKb",PARAMETER_TYPE,GKbIndex,0.003,0.0018,0.003,"mS/uF"}
};
typedef struct parameters_str { double  GKb;} PARAMETERS;
void OHaraRudy_IKbFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void OHaraRudy_IKbAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case GKbIndex:
            *value = parms->GKb; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case GKbIndex:
            parms->GKb = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO OHaraRudy_IKbInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 1;
   info.varInfo = varInfo;
   info.func = OHaraRudy_IKbFunc;
   info.access = OHaraRudy_IKbAccess;
   return info;
}
#ifdef doEnd
#define ENDCODE() {\
   if (derived->dState != 0) \
   {\
   double  *dState = derived->dState+pOffset;\
   }\
}
#else
#define ENDCODE()
#endif
