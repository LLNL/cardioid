/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ PCabIndex, nVar};
static const char *compName = "ICab";
static VARINFO varInfo[] =
{
   {"PCab",PARAMETER_TYPE,PCabIndex,2.5e-08,2.5e-08,2.5e-08,"cm/s"}
};
typedef struct parameters_str { double  PCab;} PARAMETERS;
void OHaraRudy_ICabFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void OHaraRudy_ICabAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case PCabIndex:
            *value = parms->PCab; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case PCabIndex:
            parms->PCab = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO OHaraRudy_ICabInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 1;
   info.varInfo = varInfo;
   info.func = OHaraRudy_ICabFunc;
   info.access = OHaraRudy_ICabAccess;
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
