#include <assert.h>
enum enumIndex{ PCabIndex, nVar};
static VARINFO varInfo[] =
{
   {"PCab",PARAMETER_TYPE,PCabIndex,2.5e-08,2.5e-08,2.5e-08,"cm/s"}
};
typedef struct parameters_str { double  PCab;} PARAMETERS;
void OHaraRudy_ICabAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case PCabIndex:
            *value = parms->PCab; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case PCabIndex:
            parms->PCab = *value;
            break;
            assert(0); 
      }
   }
}
