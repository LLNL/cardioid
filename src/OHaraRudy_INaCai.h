/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ GNaCaiIndex, nVar};
static const char *compName = "INaCai";
static VARINFO varInfo[] =
{
   {"GNaCai",PARAMETER_TYPE,GNaCaiIndex,0.00064,0.000704,0.000896,"uA/uF"}
};
typedef struct parameters_str { double  GNaCai;} PARAMETERS;
void OHaraRudy_INaCaiFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void OHaraRudy_INaCaiAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case GNaCaiIndex:
            *value = parms->GNaCai; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case GNaCaiIndex:
            parms->GNaCai = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO OHaraRudy_INaCaiInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 1;
   info.varInfo = varInfo;
   info.func = OHaraRudy_INaCaiFunc;
   info.access = OHaraRudy_INaCaiAccess;
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
