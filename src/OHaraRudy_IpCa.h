#include <assert.h>
enum enumIndex{ GpCaIndex, nVar};
static VARINFO varInfo[] =
{
   {"GpCa",PARAMETER_TYPE,GpCaIndex,0.0005,0.0005,0.0005,"mS/uF"}
};
typedef struct parameters_str { double  GpCa;} PARAMETERS;
void OHaraRudy_IpCaAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case GpCaIndex:
            *value = parms->GpCa; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case GpCaIndex:
            parms->GpCa = *value;
            break;
            assert(0); 
      }
   }
}
