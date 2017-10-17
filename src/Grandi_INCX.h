/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ AFIndex, GNCXIndex, nVar};
static const char *compName = "INCX";
static VARINFO varInfo[] =
{
   {"AF",PARAMETER_TYPE,AFIndex,0,0,1,1,"1"},
   {"GNCX",PARAMETER_TYPE,GNCXIndex,3.15,3.15,3.15,3.15,"uA/uF"}
};
typedef struct parameters_str { double  AF, GNCX;} PARAMETERS;
void Grandi_INCXFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void Grandi_INCXAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case AFIndex:
            *value = parms->AF; 
            break;
         case GNCXIndex:
            *value = parms->GNCX; 
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
         case GNCXIndex:
            parms->GNCX = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO Grandi_INCXInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 2;
   info.varInfo = varInfo;
   info.func = Grandi_INCXFunc;
   info.access = Grandi_INCXAccess;
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
