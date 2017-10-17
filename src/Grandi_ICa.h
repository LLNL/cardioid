/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ AFIndex, pCaIndex, pNaIndex, pKIndex, dIndex, fIndex, fcaBjIndex, fcaBslIndex, nVar};
static const char *compName = "ICa";
static VARINFO varInfo[] =
{
   {"AF",PARAMETER_TYPE,AFIndex,0,0,1,1,"1"},
   {"pCa",PARAMETER_TYPE,pCaIndex,0.00027,0.00027,0.00027,0.00027,"cm/s"},
   {"pNa",PARAMETER_TYPE,pNaIndex,7.5e-09,7.5e-09,7.5e-09,7.5e-09,"cm/s"},
   {"pK",PARAMETER_TYPE,pKIndex,1.35e-07,1.35e-07,1.35e-07,1.35e-07,"cm/s"},
   {"d",PSTATE_TYPE,dIndex,0,0,0,0,"1"},
   {"f",PSTATE_TYPE,fIndex,1,1,1,1,"1"},
   {"fcaBj",PSTATE_TYPE,fcaBjIndex,0.025,0.025,0.025,0.025,"1"},
   {"fcaBsl",PSTATE_TYPE,fcaBslIndex,0.015,0.015,0.015,0.015,"1"}
};
typedef struct parameters_str { double  AF, pCa, pNa, pK;} PARAMETERS;
typedef struct pstate_str { double  d, f, fcaBj, fcaBsl;} PSTATE;
void Grandi_ICaFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void Grandi_ICaAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PSTATE *state = (PSTATE *)statePtr;
   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case AFIndex:
            *value = parms->AF; 
            break;
         case pCaIndex:
            *value = parms->pCa; 
            break;
         case pNaIndex:
            *value = parms->pNa; 
            break;
         case pKIndex:
            *value = parms->pK; 
            break;
         case dIndex:
            *value = state->d; 
            break;
         case fIndex:
            *value = state->f; 
            break;
         case fcaBjIndex:
            *value = state->fcaBj; 
            break;
         case fcaBslIndex:
            *value = state->fcaBsl; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case AFIndex:
            parms->AF = *value;
            break;
         case pCaIndex:
            parms->pCa = *value;
            break;
         case pNaIndex:
            parms->pNa = *value;
            break;
         case pKIndex:
            parms->pK = *value;
            break;
         case dIndex:
            state->d = *value;
            break;
         case fIndex:
            state->f = *value;
            break;
         case fcaBjIndex:
            state->fcaBj = *value;
            break;
         case fcaBslIndex:
            state->fcaBsl = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO Grandi_ICaInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.pStateSize = sizeof(PSTATE);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 8;
   info.varInfo = varInfo;
   info.func = Grandi_ICaFunc;
   info.access = Grandi_ICaAccess;
   return info;
}
#ifdef doEnd
#define ENDCODE() {\
   if (derived->dState != 0) \
   {\
   double  *dState = derived->dState+pOffset;\
   dState[0]=dd; \
   dState[1]=df; \
   dState[2]=dfcaBj; \
   dState[3]=dfcaBsl; \
   }\
}
#else
#define ENDCODE()
#endif
