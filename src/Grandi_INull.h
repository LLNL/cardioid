/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ GrIIndex, nVar};
static const char *compName = "INull";
static VARINFO varInfo[] =
{
  {"GrI",PARAMETER_TYPE,GrIIndex,0,0,0,0,"mS/uF"}
};
typedef struct parameters_str { double  GrI;} PARAMETERS;
void Grandi_INullFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void Grandi_INullAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case GrIIndex:
            *value = parms->GrI; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case GrIIndex:
            parms->GrI = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO Grandi_INullInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 1;
   info.varInfo = varInfo;
   info.func = Grandi_INullFunc;
   info.access = Grandi_INullAccess;
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
