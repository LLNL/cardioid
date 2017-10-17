/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ PCaLIndex, PCaNaIndex, PCaKIndex, dIndex, fFastIndex, fSlowIndex, fCaFastIndex, fCaSlowIndex, jCaIndex, nIndex, fCaMKFastIndex, fCaCaMKFastIndex, nVar};
static const char *compName = "ICa";
static VARINFO varInfo[] =
{
   {"PCaL",PARAMETER_TYPE,PCaLIndex,0.0001,0.00012,0.00025,"cm/s"},
   {"PCaNa",PARAMETER_TYPE,PCaNaIndex,1.25e-07,1.5e-07,3.125e-07,"cm/s"},
   {"PCaK",PARAMETER_TYPE,PCaKIndex,3.574e-08,4.2888e-08,8.935e-08,"cm/s"},
   {"d",PSTATE_TYPE,dIndex,2.43015e-09,2.43015e-09,2.43015e-09,"1"},
   {"fFast",PSTATE_TYPE,fFastIndex,1,1,1,"1"},
   {"fSlow",PSTATE_TYPE,fSlowIndex,0.910671,0.910671,0.910671,"1"},
   {"fCaFast",PSTATE_TYPE,fCaFastIndex,1,1,1,"1"},
   {"fCaSlow",PSTATE_TYPE,fCaSlowIndex,0.99982,0.99982,0.99982,"1"},
   {"jCa",PSTATE_TYPE,jCaIndex,0.999977,0.999977,0.999977,"1"},
   {"n",PSTATE_TYPE,nIndex,0.00267171,0.00267171,0.00267171,"1"},
   {"fCaMKFast",PSTATE_TYPE,fCaMKFastIndex,1,1,1,"1"},
   {"fCaCaMKFast",PSTATE_TYPE,fCaCaMKFastIndex,1,1,1,"1"}
};
typedef struct parameters_str { double  PCaL, PCaNa, PCaK;} PARAMETERS;
typedef struct pstate_str { double  d, fFast, fSlow, fCaFast, fCaSlow, jCa, n, fCaMKFast, fCaCaMKFast;} PSTATE;
void OHaraRudy_ICaFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void OHaraRudy_ICaAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PSTATE *state = (PSTATE *)statePtr;
   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case PCaLIndex:
            *value = parms->PCaL; 
            break;
         case PCaNaIndex:
            *value = parms->PCaNa; 
            break;
         case PCaKIndex:
            *value = parms->PCaK; 
            break;
         case dIndex:
            *value = state->d; 
            break;
         case fFastIndex:
            *value = state->fFast; 
            break;
         case fSlowIndex:
            *value = state->fSlow; 
            break;
         case fCaFastIndex:
            *value = state->fCaFast; 
            break;
         case fCaSlowIndex:
            *value = state->fCaSlow; 
            break;
         case jCaIndex:
            *value = state->jCa; 
            break;
         case nIndex:
            *value = state->n; 
            break;
         case fCaMKFastIndex:
            *value = state->fCaMKFast; 
            break;
         case fCaCaMKFastIndex:
            *value = state->fCaCaMKFast; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case PCaLIndex:
            parms->PCaL = *value;
            break;
         case PCaNaIndex:
            parms->PCaNa = *value;
            break;
         case PCaKIndex:
            parms->PCaK = *value;
            break;
         case dIndex:
            state->d = *value;
            break;
         case fFastIndex:
            state->fFast = *value;
            break;
         case fSlowIndex:
            state->fSlow = *value;
            break;
         case fCaFastIndex:
            state->fCaFast = *value;
            break;
         case fCaSlowIndex:
            state->fCaSlow = *value;
            break;
         case jCaIndex:
            state->jCa = *value;
            break;
         case nIndex:
            state->n = *value;
            break;
         case fCaMKFastIndex:
            state->fCaMKFast = *value;
            break;
         case fCaCaMKFastIndex:
            state->fCaCaMKFast = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO OHaraRudy_ICaInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.pStateSize = sizeof(PSTATE);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 12;
   info.varInfo = varInfo;
   info.func = OHaraRudy_ICaFunc;
   info.access = OHaraRudy_ICaAccess;
   return info;
}
#ifdef doEnd
#define ENDCODE() {\
   if (derived->dState != 0) \
   {\
   double  *dState = derived->dState+pOffset;\
   dState[0]=dd; \
   dState[1]=dfFast; \
   dState[2]=dfSlow; \
   dState[3]=dfCaFast; \
   dState[4]=dfCaSlow; \
   dState[5]=djCa; \
   dState[6]=dn; \
   dState[7]=dfCaMKFast; \
   dState[8]=dfCaCaMKFast; \
   }\
}
#else
#define ENDCODE()
#endif
