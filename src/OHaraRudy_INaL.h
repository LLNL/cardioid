#include <assert.h>
enum enumIndex{ GNaLIndex, mLIndex, hLIndex, hLCaMKIndex, nVar};
static VARINFO varInfo[] =
{
   {"GNaL",PARAMETER_TYPE,GNaLIndex,0.0075,0.0045,0.0075,"mS/uF"},
   {"mL",PSTATE_TYPE,mLIndex,0.000194015,0.000194015,0.000194015,"1"},
   {"hL",PSTATE_TYPE,hLIndex,0.496116,0.496116,0.496116,"1"},
   {"hLCaMK",PSTATE_TYPE,hLCaMKIndex,0.265885,0.265885,0.265885,"1"}
};
typedef struct parameters_str { double  GNaL;} PARAMETERS;
typedef struct pstate_str { double  mL, hL, hLCaMK;} PSTATE;
void OHaraRudy_INaLAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PSTATE *state = (PSTATE *)statePtr;
   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case GNaLIndex:
            *value = parms->GNaL; 
            break;
         case mLIndex:
            *value = state->mL; 
            break;
         case hLIndex:
            *value = state->hL; 
            break;
         case hLCaMKIndex:
            *value = state->hLCaMK; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case GNaLIndex:
            parms->GNaL = *value;
            break;
         case mLIndex:
            state->mL = *value;
            break;
         case hLIndex:
            state->hL = *value;
            break;
         case hLCaMKIndex:
            state->hLCaMK = *value;
            break;
            assert(0); 
      }
   }
}
