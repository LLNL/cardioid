#include <assert.h>
enum enumIndex{ GNaCassIndex, nVar};
static VARINFO varInfo[] =
{
   {"GNaCass",PARAMETER_TYPE,GNaCassIndex,0.00016,0.000176,0.000224,"uA/uF"}
};
typedef struct parameters_str { double  GNaCass;} PARAMETERS;
void OHaraRudy_INaCassAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case GNaCassIndex:
            *value = parms->GNaCass; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case GNaCassIndex:
            parms->GNaCass = *value;
            break;
            assert(0); 
      }
   }
}
