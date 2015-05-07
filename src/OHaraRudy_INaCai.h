#include <assert.h>
enum enumIndex{ GNaCaiIndex, nVar};
static VARINFO varInfo[] =
{
   {"GNaCai",PARAMETER_TYPE,GNaCaiIndex,0.00064,0.000704,0.000896,"uA/uF"}
};
typedef struct parameters_str { double  GNaCai;} PARAMETERS;
void OHaraRudy_INaCaiAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case GNaCaiIndex:
            *value = parms->GNaCai; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case GNaCaiIndex:
            parms->GNaCai = *value;
            break;
            assert(0); 
      }
   }
}
