/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ PNabIndex, nVar};
static const char *compName = "INab";
static VARINFO varInfo[] =
{
   {"PNab",PARAMETER_TYPE,PNabIndex,3.75e-10,3.75e-10,3.75e-10,"cm/s"}
};
typedef struct parameters_str { double  PNab;} PARAMETERS;
void OHaraRudy_INabFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void OHaraRudy_INabAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case PNabIndex:
            *value = parms->PNab; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case PNabIndex:
            parms->PNab = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO OHaraRudy_INabInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 1;
   info.varInfo = varInfo;
   info.func = OHaraRudy_INabFunc;
   info.access = OHaraRudy_INabAccess;
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
