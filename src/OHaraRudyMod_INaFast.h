/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ GNaFastIndex, mIndex, hFastIndex, hSlowIndex, jIndex, hCaMKSlowIndex, jCaMKIndex, nVar};
static const char *compName = "INaFast";
static VARINFO varInfo[] =
{
   {"GNaFast",PARAMETER_TYPE,GNaFastIndex,75,75,75,"mS/uF"},
   {"m",PSTATE_TYPE,mIndex,0.0074621,0.0074621,0.0074621,"1"},
   {"hFast",PSTATE_TYPE,hFastIndex,0.692591,0.692591,0.692591,"1"},
   {"hSlow",PSTATE_TYPE,hSlowIndex,0.692574,0.692574,0.692574,"1"},
   {"j",PSTATE_TYPE,jIndex,0.692477,0.692477,0.692477,"1"},
   {"hCaMKSlow",PSTATE_TYPE,hCaMKSlowIndex,0.448501,0.448501,0.448501,"1"},
   {"jCaMK",PSTATE_TYPE,jCaMKIndex,0.692413,0.692413,0.692413,"1"}
};
typedef struct parameters_str { double  GNaFast;} PARAMETERS;
typedef struct pstate_str { double  m, hFast, hSlow, j, hCaMKSlow, jCaMK;} PSTATE;
void OHaraRudyMod_INaFastFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void OHaraRudyMod_INaFastAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PSTATE *state = (PSTATE *)statePtr;
   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case GNaFastIndex:
            *value = parms->GNaFast; 
            break;
         case mIndex:
            *value = state->m; 
            break;
         case hFastIndex:
            *value = state->hFast; 
            break;
         case hSlowIndex:
            *value = state->hSlow; 
            break;
         case jIndex:
            *value = state->j; 
            break;
         case hCaMKSlowIndex:
            *value = state->hCaMKSlow; 
            break;
         case jCaMKIndex:
            *value = state->jCaMK; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case GNaFastIndex:
            parms->GNaFast = *value;
            break;
         case mIndex:
            state->m = *value;
            break;
         case hFastIndex:
            state->hFast = *value;
            break;
         case hSlowIndex:
            state->hSlow = *value;
            break;
         case jIndex:
            state->j = *value;
            break;
         case hCaMKSlowIndex:
            state->hCaMKSlow = *value;
            break;
         case jCaMKIndex:
            state->jCaMK = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO OHaraRudyMod_INaFastInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.pStateSize = sizeof(PSTATE);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 7;
   info.varInfo = varInfo;
   info.func = OHaraRudyMod_INaFastFunc;
   info.access = OHaraRudyMod_INaFastAccess;
   return info;
}
#ifdef doEnd
#define ENDCODE() {\
   if (derived->dState != 0) \
   {\
   double  *dState = derived->dState+pOffset;\
   dState[0]=dm; \
   dState[1]=dhFast; \
   dState[2]=dhSlow; \
   dState[3]=dj; \
   dState[4]=dhCaMKSlow; \
   dState[5]=djCaMK; \
   }\
}
#else
#define ENDCODE()
#endif
