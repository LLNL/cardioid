/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ GClCaIndex, nVar};
static const char *compName = "IClCa";
static VARINFO varInfo[] =
{
   {"GClCa",PARAMETER_TYPE,GClCaIndex,0.0548,0.0548,0.0548,0.0548,"mS/uF"}
};
typedef struct parameters_str { double  GClCa;} PARAMETERS;
void Grandi_IClCaFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void Grandi_IClCaAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case GClCaIndex:
            *value = parms->GClCa; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case GClCaIndex:
            parms->GClCa = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO Grandi_IClCaInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 1;
   info.varInfo = varInfo;
   info.func = Grandi_IClCaFunc;
   info.access = Grandi_IClCaAccess;
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
