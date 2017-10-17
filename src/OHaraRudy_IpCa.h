/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ GpCaIndex, nVar};
static const char *compName = "IpCa";
static VARINFO varInfo[] =
{
   {"GpCa",PARAMETER_TYPE,GpCaIndex,0.0005,0.0005,0.0005,"mS/uF"}
};
typedef struct parameters_str { double  GpCa;} PARAMETERS;
void OHaraRudy_IpCaFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void OHaraRudy_IpCaAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case GpCaIndex:
            *value = parms->GpCa; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case GpCaIndex:
            parms->GpCa = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO OHaraRudy_IpCaInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 1;
   info.varInfo = varInfo;
   info.func = OHaraRudy_IpCaFunc;
   info.access = OHaraRudy_IpCaAccess;
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
