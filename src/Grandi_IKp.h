/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ GKpIndex, nVar};
static const char *compName = "IKp";
static VARINFO varInfo[] =
{
   {"GKp",PARAMETER_TYPE,GKpIndex,0.002,0.002,0.002,0.002,"mS/uF"}
};
typedef struct parameters_str { double  GKp;} PARAMETERS;
void Grandi_IKpFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void Grandi_IKpAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case GKpIndex:
            *value = parms->GKp; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case GKpIndex:
            parms->GKp = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO Grandi_IKpInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 1;
   info.varInfo = varInfo;
   info.func = Grandi_IKpFunc;
   info.access = Grandi_IKpAccess;
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
