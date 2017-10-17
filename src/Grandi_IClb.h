/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ GClbIndex, nVar};
static const char *compName = "IClb";
static VARINFO varInfo[] =
{
   {"GClb",PARAMETER_TYPE,GClbIndex,0.009,0.009,0.009,0.009,"mS/uF"}
};
typedef struct parameters_str { double  GClb;} PARAMETERS;
void Grandi_IClbFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void Grandi_IClbAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case GClbIndex:
            *value = parms->GClb; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case GClbIndex:
            parms->GClb = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO Grandi_IClbInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 1;
   info.varInfo = varInfo;
   info.func = Grandi_IClbFunc;
   info.access = Grandi_IClbAccess;
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
