/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ VmIndex, dVmIndex, iStimIndex, nVar};
static const char *compName = "Voltage";
static VARINFO varInfo[] =
{
   {"Vm",PSTATE_TYPE,VmIndex,-87.84,-87.84,-87.84,"mV"},
   {"dVm",PSTATE_TYPE,dVmIndex,0,0,0,"mV/ms"},
   {"iStim",PSTATE_TYPE,iStimIndex,0,0,0,"mV/ms"}
};
typedef struct pstate_str { double  Vm, dVm, iStim;} PSTATE;
void OHaraRudy_VoltageFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void OHaraRudy_VoltageAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PSTATE *state = (PSTATE *)statePtr;
   if (type == READ)
   {
      switch (index)
      {
         case VmIndex:
            *value = state->Vm; 
            break;
         case dVmIndex:
            *value = state->dVm; 
            break;
         case iStimIndex:
            *value = state->iStim; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case VmIndex:
            state->Vm = *value;
            break;
         case dVmIndex:
            state->dVm = *value;
            break;
         case iStimIndex:
            state->iStim = *value;
            break;
            assert(0); 
      }
   }
}
COMPONENTINFO OHaraRudy_VoltageInit()
{
   COMPONENTINFO info;
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.pStateSize = sizeof(PSTATE);
   info.nVar = 3;
   info.varInfo = varInfo;
   info.func = OHaraRudy_VoltageFunc;
   info.access = OHaraRudy_VoltageAccess;
   return info;
}
#ifdef doEnd
#define ENDCODE() {\
   if (derived->dState != 0) \
   {\
   double  *dState = derived->dState+pOffset;\
   dState[0]=dVm; \
   dState[1]=ddVm; \
   dState[2]=diStim; \
   }\
}
#else
#define ENDCODE()
#endif
