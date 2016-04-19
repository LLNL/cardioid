/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ IbarSLCaPIndex, nVar};
static const char *compName = "IpCa";
static VARINFO varInfo[] =
{
  {"IbarSLCaP",PARAMETER_TYPE,IbarSLCaPIndex,0.0471,0.0471,0.0471,0.0471,"mS/uF"}
};
typedef struct parameters_str { double  IbarSLCaP;} PARAMETERS;
void Grandi_IpCaFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void Grandi_IpCaAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case IbarSLCaPIndex:
            *value = parms->IbarSLCaP; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case IbarSLCaPIndex:
            parms->IbarSLCaP = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO Grandi_IpCaInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 1;
   info.varInfo = varInfo;
   info.func = Grandi_IpCaFunc;
   info.access = Grandi_IpCaAccess;
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
