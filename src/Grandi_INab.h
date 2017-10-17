/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ GNabIndex, nVar};
static const char *compName = "INab";
static VARINFO varInfo[] =
{
   {"GNab",PARAMETER_TYPE,GNabIndex,0.000597,0.000597,0.000597,0.000597,"mS/uF"}
};
typedef struct parameters_str { double  GNab;} PARAMETERS;
void Grandi_INabFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void Grandi_INabAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case GNabIndex:
            *value = parms->GNab; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case GNabIndex:
            parms->GNab = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO Grandi_INabInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 1;
   info.varInfo = varInfo;
   info.func = Grandi_INabFunc;
   info.access = Grandi_INabAccess;
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
