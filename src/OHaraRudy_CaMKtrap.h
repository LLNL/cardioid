#include <assert.h>
enum enumIndex{ CaMKtrapIndex, nVar};
static VARINFO varInfo[] =
{
   {"CaMKtrap",PSTATE_TYPE,CaMKtrapIndex,0.0124065,0.0124065,0.0124065,"1"}
};
typedef struct pstate_str { double  CaMKtrap;} PSTATE;
void OHaraRudy_CaMKtrapAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PSTATE *state = (PSTATE *)statePtr;
   if (type == READ)
   {
      switch (index)
      {
         case CaMKtrapIndex:
            *value = state->CaMKtrap; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case CaMKtrapIndex:
            state->CaMKtrap = *value;
            break;
            assert(0); 
      }
   }
}
