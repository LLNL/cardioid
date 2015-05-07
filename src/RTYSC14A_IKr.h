#include <assert.h>
enum enumIndex{ GKrIndex, C3Index, C2Index, C1Index, OIndex, IIndex, nVar};
static VARINFO varInfo[] =
{
   {"GKr",PARAMETER_TYPE,GKrIndex,0.046,0.0598,0.0368,"mS/uF"},
   {"C3",PSTATE_TYPE,C3Index,1,1,1,"1"},
   {"C2",PSTATE_TYPE,C2Index,0,0,0,"1"},
   {"C1",PSTATE_TYPE,C1Index,0,0,0,"1"},
   {"O",PSTATE_TYPE,OIndex,0,0,0,"1"},
   {"I",PSTATE_TYPE,IIndex,0,0,0,"1"}
};
typedef struct parameters_str { double  GKr;} PARAMETERS;
typedef struct pstate_str { double  C3, C2, C1, O, I;} PSTATE;
void RTYSC14A_IKrAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
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
         case C3Index:
            *value = state->C3; 
            break;
         case C2Index:
            *value = state->C2; 
            break;
         case C1Index:
            *value = state->C1; 
            break;
         case OIndex:
            *value = state->O; 
            break;
         case IIndex:
            *value = state->I; 
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
         case C3Index:
            state->C3 = *value;
            break;
         case C2Index:
            state->C2 = *value;
            break;
         case C1Index:
            state->C1 = *value;
            break;
         case OIndex:
            state->O = *value;
            break;
         case IIndex:
            state->I = *value;
            break;
            assert(0); 
      }
   }
}
