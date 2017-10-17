/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ GCabIndex, nVar};
static const char *compName = "ICab";
static VARINFO varInfo[] =
{
   {"GCab",PARAMETER_TYPE,GCabIndex,0.00060643,0.00060643,0.00060643,0.00060643,"mS/uF"}
};
typedef struct parameters_str { double  GCab;} PARAMETERS;
void Grandi_ICabFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void Grandi_ICabAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case GCabIndex:
            *value = parms->GCab; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case GCabIndex:
            parms->GCab = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO Grandi_ICabInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 1;
   info.varInfo = varInfo;
   info.func = Grandi_ICabFunc;
   info.access = Grandi_ICabAccess;
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
