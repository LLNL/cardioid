/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */
#include <assert.h>
#include <string.h>
enum enumIndex{ CMDNIndex, NaiIndex, NassIndex, KiIndex, KssIndex, CaiIndex, CassIndex, CansrIndex, CajsrIndex, nVar};
static const char *compName = "Concent";
static VARINFO varInfo[] =
{
   {"CMDN",PARAMETER_TYPE,CMDNIndex,0.05,0.065,0.05,"mM"},
   {"Nai",PSTATE_TYPE,NaiIndex,7.23,7.23,7.23,"mM"},
   {"Nass",PSTATE_TYPE,NassIndex,7.23,7.23,7.23,"mM"},
   {"Ki",PSTATE_TYPE,KiIndex,143.79,143.79,143.79,"mM"},
   {"Kss",PSTATE_TYPE,KssIndex,143.79,143.79,143.79,"mM"},
   {"Cai",PSTATE_TYPE,CaiIndex,8.54e-05,8.54e-05,8.54e-05,"mM"},
   {"Cass",PSTATE_TYPE,CassIndex,8.43e-05,8.43e-05,8.43e-05,"mM"},
   {"Cansr",PSTATE_TYPE,CansrIndex,1.61,1.61,1.61,"mM"},
   {"Cajsr",PSTATE_TYPE,CajsrIndex,1.56,1.56,1.56,"mM"}
};
typedef struct parameters_str { double  CMDN;} PARAMETERS;
typedef struct pstate_str { double  Nai, Nass, Ki, Kss, Cai, Cass, Cansr, Cajsr;} PSTATE;
void OHaraRudy_ConcentFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );
void OHaraRudy_ConcentAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)
{

   PSTATE *state = (PSTATE *)statePtr;
   PARAMETERS *parms = (PARAMETERS *)parmsPtr;
   if (type == READ)
   {
      switch (index)
      {
         case CMDNIndex:
            *value = parms->CMDN; 
            break;
         case NaiIndex:
            *value = state->Nai; 
            break;
         case NassIndex:
            *value = state->Nass; 
            break;
         case KiIndex:
            *value = state->Ki; 
            break;
         case KssIndex:
            *value = state->Kss; 
            break;
         case CaiIndex:
            *value = state->Cai; 
            break;
         case CassIndex:
            *value = state->Cass; 
            break;
         case CansrIndex:
            *value = state->Cansr; 
            break;
         case CajsrIndex:
            *value = state->Cajsr; 
            break;
         default:
            assert(0); 
      }
   }
   if (type == WRITE)
   {
      switch (index)
      {
         case CMDNIndex:
            parms->CMDN = *value;
            break;
         case NaiIndex:
            state->Nai = *value;
            break;
         case NassIndex:
            state->Nass = *value;
            break;
         case KiIndex:
            state->Ki = *value;
            break;
         case KssIndex:
            state->Kss = *value;
            break;
         case CaiIndex:
            state->Cai = *value;
            break;
         case CassIndex:
            state->Cass = *value;
            break;
         case CansrIndex:
            state->Cansr = *value;
            break;
         case CajsrIndex:
            state->Cajsr = *value;
            break;
            assert(0); 
      }
   }
}
void   OHaraRudy_ConcentConstants();
COMPONENTINFO OHaraRudy_ConcentInit()
{
   COMPONENTINFO info;
   OHaraRudy_ConcentConstants();
   if (FRT  < 0) FRT = F/(R*T);
   info.compName = strdup(compName);
   info.pStateSize = sizeof(PSTATE);
   info.parmsSize = sizeof(PARAMETERS);
   info.nVar = 9;
   info.varInfo = varInfo;
   info.func = OHaraRudy_ConcentFunc;
   info.access = OHaraRudy_ConcentAccess;
   return info;
}
#ifdef doEnd
#define ENDCODE() {\
   if (derived->dState != 0) \
   {\
   double  *dState = derived->dState+pOffset;\
   dState[0]=dNai; \
   dState[1]=dNass; \
   dState[2]=dKi; \
   dState[3]=dKss; \
   dState[4]=dCai; \
   dState[5]=dCass; \
   dState[6]=dCansr; \
   dState[7]=dCajsr; \
   }\
}
#else
#define ENDCODE()
#endif
