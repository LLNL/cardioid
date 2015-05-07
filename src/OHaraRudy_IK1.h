#include <assert.h>
enum enumIndex{ GK1Index, XK1Index, nVar};
static VARINFO varInfo[] =
{
   {"GK1",PARAMETER_TYPE,GK1Index,0.1908,0.22896,0.24804,"mS/uF"},
   {"XK1",PSTATE_TYPE,XK1Index,0.996801,0.996801,0.996801,"1"}
};
typedef struct parameters_str { double  GK1;} PARAMETERS;
typedef struct pstate_str { double  XK1;} PSTATE;
void OHaraRudy_IK1Access(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PSTATE *state = (PSTATE *)statePtr;
   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case GK1Index:
            *value = parms->GK1; 
            break;
         case XK1Index:
            *value = state->XK1; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case GK1Index:
            parms->GK1 = *value;
            break;
         case XK1Index:
            state->XK1 = *value;
            break;
            assert(0); 
      }
   }
}
