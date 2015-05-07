#include <assert.h>
enum enumIndex{ PNabIndex, nVar};
static VARINFO varInfo[] =
{
   {"PNab",PARAMETER_TYPE,PNabIndex,3.75e-10,3.75e-10,3.75e-10,"cm/s"}
};
typedef struct parameters_str { double  PNab;} PARAMETERS;
void OHaraRudy_INabAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case PNabIndex:
            *value = parms->PNab; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case PNabIndex:
            parms->PNab = *value;
            break;
            assert(0); 
      }
   }
}
