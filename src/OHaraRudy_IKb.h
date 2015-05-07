#include <assert.h>
enum enumIndex{ GKbIndex, nVar};
static VARINFO varInfo[] =
{
   {"GKb",PARAMETER_TYPE,GKbIndex,0.003,0.0018,0.003,"mS/uF"}
};
typedef struct parameters_str { double  GKb;} PARAMETERS;
void OHaraRudy_IKbAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case GKbIndex:
            *value = parms->GKb; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case GKbIndex:
            parms->GKb = *value;
            break;
            assert(0); 
      }
   }
}
