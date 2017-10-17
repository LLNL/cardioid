/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ PNaKIndex, nVar};
static const char *compName = "INaK";
static VARINFO varInfo[] =
{
   {"PNaK",PARAMETER_TYPE,PNaKIndex,30,27,21,"mV/mM"}
};
typedef struct parameters_str { double  PNaK;} PARAMETERS;
void OHaraRudy_INaKFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void OHaraRudy_INaKAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case PNaKIndex:
            *value = parms->PNaK; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case PNaKIndex:
            parms->PNaK = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO OHaraRudy_INaKInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 1;
   info.varInfo = varInfo;
   info.func = OHaraRudy_INaKFunc;
   info.access = OHaraRudy_INaKAccess;
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
