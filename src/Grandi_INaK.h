/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ IbarNaKIndex, nVar};
static const char *compName = "INaK";
static VARINFO varInfo[] =
{
  {"IbarNaK",PARAMETER_TYPE,IbarNaKIndex,1.26,1.26,1.26,1.26,"mV/mM"}
};
typedef struct parameters_str { double  IbarNaK;} PARAMETERS;
void Grandi_INaKFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void Grandi_INaKAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case IbarNaKIndex:
            *value = parms->IbarNaK; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case IbarNaKIndex:
            parms->IbarNaK = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO Grandi_INaKInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 1;
   info.varInfo = varInfo;
   info.func = Grandi_INaKFunc;
   info.access = Grandi_INaKAccess;
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
