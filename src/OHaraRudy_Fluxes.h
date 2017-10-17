/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ alphaJrelNPIndex, betaJrelNPIndex, alphaJrelCaMKIndex, betaJrelCaMKIndex, cJupIndex, JrelNPIndex, JrelCaMKIndex, nVar};
static const char *compName = "Fluxes";
static VARINFO varInfo[] =
{
   {"alphaJrelNP",PARAMETER_TYPE,alphaJrelNPIndex,2.375,2.375,4.0375,"mM/mV"},
   {"betaJrelNP",PARAMETER_TYPE,betaJrelNPIndex,4.75,4.75,4.75,"mM/mS"},
   {"alphaJrelCaMK",PARAMETER_TYPE,alphaJrelCaMKIndex,2.96875,2.96875,5.046875,"mM/mV"},
   {"betaJrelCaMK",PARAMETER_TYPE,betaJrelCaMKIndex,5.9375,5.9375,5.9375,"mM/mS"},
   {"cJup",PARAMETER_TYPE,cJupIndex,0.004375,0.0056875,0.004375,"mM/ms"},
   {"JrelNP",PSTATE_TYPE,JrelNPIndex,2.53943e-05,2.53943e-05,2.53943e-05,"mM/ms"},
   {"JrelCaMK",PSTATE_TYPE,JrelCaMKIndex,3.17262e-07,3.17262e-07,3.17262e-07,"mM/ms"}
};
typedef struct parameters_str { double  alphaJrelNP, betaJrelNP, alphaJrelCaMK, betaJrelCaMK, cJup;} PARAMETERS;
typedef struct pstate_str { double  JrelNP, JrelCaMK;} PSTATE;
void OHaraRudy_FluxesFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void OHaraRudy_FluxesAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PSTATE *state = (PSTATE *)statePtr;
   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case alphaJrelNPIndex:
            *value = parms->alphaJrelNP; 
            break;
         case betaJrelNPIndex:
            *value = parms->betaJrelNP; 
            break;
         case alphaJrelCaMKIndex:
            *value = parms->alphaJrelCaMK; 
            break;
         case betaJrelCaMKIndex:
            *value = parms->betaJrelCaMK; 
            break;
         case cJupIndex:
            *value = parms->cJup; 
            break;
         case JrelNPIndex:
            *value = state->JrelNP; 
            break;
         case JrelCaMKIndex:
            *value = state->JrelCaMK; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case alphaJrelNPIndex:
            parms->alphaJrelNP = *value;
            break;
         case betaJrelNPIndex:
            parms->betaJrelNP = *value;
            break;
         case alphaJrelCaMKIndex:
            parms->alphaJrelCaMK = *value;
            break;
         case betaJrelCaMKIndex:
            parms->betaJrelCaMK = *value;
            break;
         case cJupIndex:
            parms->cJup = *value;
            break;
         case JrelNPIndex:
            state->JrelNP = *value;
            break;
         case JrelCaMKIndex:
            state->JrelCaMK = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO OHaraRudy_FluxesInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.pStateSize = sizeof(PSTATE);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 7;
   info.varInfo = varInfo;
   info.func = OHaraRudy_FluxesFunc;
   info.access = OHaraRudy_FluxesAccess;
   return info;
}
#ifdef doEnd
#define ENDCODE() {\
   if (derived->dState != 0) \
   {\
   double  *dState = derived->dState+pOffset;\
   dState[0]=dJrelNP; \
   dState[1]=dJrelCaMK; \
   }\
}
#else
#define ENDCODE()
#endif
