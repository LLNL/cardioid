/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ AFIndex, GK1Index, nVar};
static const char *compName = "IK1";
static VARINFO varInfo[] =
{
   {"AF",PARAMETER_TYPE,AFIndex,0,0,1,1,"1"},
   {"GK1",PARAMETER_TYPE,GK1Index,0.0525,0.0525,0.0525,0.0525,"mS/uF"}
};
typedef struct parameters_str { double  AF, GK1;} PARAMETERS;
void Grandi_IK1Func(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void Grandi_IK1Access(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case AFIndex:
            *value = parms->AF; 
            break;
         case GK1Index:
            *value = parms->GK1; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case AFIndex:
            parms->AF = *value;
            break;
         case GK1Index:
            parms->GK1 = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO Grandi_IK1Init()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 2;
   info.varInfo = varInfo;
   info.func = Grandi_IK1Func;
   info.access = Grandi_IK1Access;
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
