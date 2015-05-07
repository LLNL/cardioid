#include <assert.h>
enum enumIndex{ GKrIndex, XrFastIndex, XrSlowIndex, nVar};
static VARINFO varInfo[] =
{
   {"GKr",PARAMETER_TYPE,GKrIndex,0.046,0.0598,0.0368,"mS/uF"},
   {"XrFast",PSTATE_TYPE,XrFastIndex,8.26608e-06,8.26608e-06,8.26608e-06,"1"},
   {"XrSlow",PSTATE_TYPE,XrSlowIndex,0.453268,0.453268,0.453268,"1"}
};
typedef struct parameters_str { double  GKr;} PARAMETERS;
typedef struct pstate_str { double  XrFast, XrSlow;} PSTATE;
void OHaraRudy_IKrAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PSTATE *state = (PSTATE *)statePtr;
   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case GKrIndex:
            *value = parms->GKr; 
            break;
         case XrFastIndex:
            *value = state->XrFast; 
            break;
         case XrSlowIndex:
            *value = state->XrSlow; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case GKrIndex:
            parms->GKr = *value;
            break;
         case XrFastIndex:
            state->XrFast = *value;
            break;
         case XrSlowIndex:
            state->XrSlow = *value;
            break;
            assert(0); 
      }
   }
}
