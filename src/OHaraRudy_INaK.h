#include <assert.h>
enum enumIndex{ PNaKIndex, nVar};
static VARINFO varInfo[] =
{
   {"PNaK",PARAMETER_TYPE,PNaKIndex,30,27,21,"mV/mM"}
};
typedef struct parameters_str { double  PNaK;} PARAMETERS;
void OHaraRudy_INaKAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case PNaKIndex:
            *value = parms->PNaK; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case PNaKIndex:
            parms->PNaK = *value;
            break;
            assert(0); 
      }
   }
}
