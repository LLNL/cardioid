/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ GNaCassIndex, nVar};
static const char *compName = "INaCass";
static VARINFO varInfo[] =
{
   {"GNaCass",PARAMETER_TYPE,GNaCassIndex,0.00016,0.000176,0.000224,"uA/uF"}
};
typedef struct parameters_str { double  GNaCass;} PARAMETERS;
void OHaraRudy_INaCassFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void OHaraRudy_INaCassAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case GNaCassIndex:
            *value = parms->GNaCass; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case GNaCassIndex:
            parms->GNaCass = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO OHaraRudy_INaCassInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 1;
   info.varInfo = varInfo;
   info.func = OHaraRudy_INaCassFunc;
   info.access = OHaraRudy_INaCassAccess;
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
