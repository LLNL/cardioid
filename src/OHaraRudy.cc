#include "OHaraRudy.hh"
#include <cmath>
#include <cassert>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <string.h>
#include "OHaraRudy.h"
#include "OHaraRudyInit.h"

using namespace std;

int OHaraRudy::parmsSize_=0; 
int OHaraRudy::stateSize_=0; 
int * OHaraRudy::privateStateOffset_=0;; 
int * OHaraRudy::privateParmsOffset_=0;; 
double * OHaraRudy::defaultStateENDO_=0;; 
double * OHaraRudy::defaultStateEPI_=0;; 
double * OHaraRudy::defaultStateM_=0;; 
double * OHaraRudy::defaultParmsENDO_=0;; 
double * OHaraRudy::defaultParmsEPI_=0;; 
double * OHaraRudy::defaultParmsM_=0;; 
COMPONENTINFO* OHaraRudy::info_; 
HandleMap  OHaraRudy::handleMap_ ;
HandleMap handleMap ;

OHaraRudy::OHaraRudy(int cellType, OHaraRudy_Parms &parms)
{

   initConsts(cellType,parms); 
   state_ = (double *)malloc(stateSize_); 
   initState(cellType);
   initParms(cellType);
}
double OHaraRudy::calc(double dt, double Vm, double iStim)
{
   DERIVED derived; 
   double Nai = ((STATE*)state_)->Nai; 
   double Ki = ((STATE*)state_)->Ki; 

   ((STATE*)state_)->Vm=Vm; 
   derived.I.stimulus = iStim; 
   reversalPotentials(Nai,Ki,&derived); 

   for (int i=1;i<=16;i++) 
      info_[i].func(cellParms_+privateParmsOffset_[i], (STATE *)state_, privateStateOffset_[ i], &derived, dt);
   info_[0].func(cellParms_+privateParmsOffset_[0], (STATE *)state_, privateStateOffset_[0], &derived, dt );

   //if (isfinite(dV)==0)MPI_Abort(MPI_COMM_WORLD,1); 
   //static double t=0; 
   //CURRENTS I=derived.I; 
   //printf("%f %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",t,I.NaCai,I.NaCass,I.NaK,I.Nab,I.Cab,I.Kb,I.pCa,I.NaFast, I.NaL, I.to, I.Kr, I.Ks,I.K1,I.CaL, I.CaNa, I.CaK); 
   //printf("IKr: %f %e\n",t,I.Kr); 
   //t += dt; 
   //exit(0); 
   return derived.dVm;  
}
double OHaraRudyDebug::calc(double dt, double Vm, double iStim)
{
   double dVm =0;
   assert(0); 
   return dVm;  
}


double OHaraRudy::defaultVoltage()
{
   return defaultVoltage_;
}

void OHaraRudy::getCheckpointInfo(vector<string>& name,
      vector<string>& unit)
{
   const HandleMap& handleMap = getHandleMap();
   for (HandleMap::const_iterator
         iter=handleMap.begin(); iter!=handleMap.end(); ++iter)
   {
      if (iter->second.checkpoint_)
      {
         name.push_back(iter->first);
         unit.push_back(iter->second.unit_);
      }
   }
}

/* the handle representation.  Returns the value undefinedName for * unrecognized varName. */
int OHaraRudy::getVarHandle(const string& varName)
{
   return getHandleMap()[varName].handle_;
}

void OHaraRudy::getValue(const vector<int>& handle,
      vector<double>& value) const
{
   for (unsigned ii=0; ii<handle.size(); ++ii)
      value[ii] = getValue(handle[ii]);
}

const string OHaraRudy::getUnit(const string& varName)
{
   return getHandleMap()[varName].unit_;
}

void OHaraRudy::initConsts(int cellType,OHaraRudy_Parms &parms)
{
   static bool initialized = false;
   if (! initialized)
   {
      typedef COMPONENTINFO (*CINIT)(); 
      int nI = parms.currentNames.size(); 
      int nComp = 3 + nI; 
      CINIT *cInit = (CINIT *)calloc(nComp,sizeof(CINIT)); 
      reversalPotentialsInit();
      int k=0; 
      cInit[k++] = OHaraRudy_ConcentInit;   // Concentration must be first and CaMK trap second. 
      cInit[k++] = OHaraRudy_CaMKtrapInit; 
      for (int i=0;i<parms.currentNames.size();i++)   // Order of currents initialization does not matter. 
      {
         string name = parms.currentModels[i] + "_" + parms.currentNames[i] ;
         if ( name == "OHaraRudy_INaCai")  cInit[k++] = OHaraRudy_INaCaiInit; 
         if ( name == "OHaraRudy_INaCass") cInit[k++] = OHaraRudy_INaCassInit; 
         if ( name == "OHaraRudy_INaK")    cInit[k++] = OHaraRudy_INaKInit; 
         if ( name == "OHaraRudy_INab")    cInit[k++] = OHaraRudy_INabInit; 
         if ( name == "OHaraRudy_ICab")    cInit[k++] = OHaraRudy_ICabInit; 
         if ( name == "OHaraRudy_IKb")     cInit[k++] = OHaraRudy_IKbInit; 
         if ( name == "OHaraRudy_IpCa")    cInit[k++] = OHaraRudy_IpCaInit; 
         if ( name == "OHaraRudy_INaFast") cInit[k++] = OHaraRudy_INaFastInit; 
         if ( name == "OHaraRudy_INaL")    cInit[k++] = OHaraRudy_INaLInit; 
         if ( name == "OHaraRudy_Ito")     cInit[k++] = OHaraRudy_ItoInit; 
         if ( name == "OHaraRudy_IKr")     cInit[k++] = OHaraRudy_IKrInit; 
         if ( name == "OHaraRudy_IKs")     cInit[k++] = OHaraRudy_IKsInit; 
         if ( name == "OHaraRudy_IK1")     cInit[k++] = OHaraRudy_IK1Init; 
         if ( name == "OHaraRudy_ICa")     cInit[k++] = OHaraRudy_ICaInit; 

         if ( name == "OHaraRudyMod_INaFast") cInit[k++] = OHaraRudyMod_INaFastInit; 
         if ( name == "RTYSC14A_IKr")      cInit[k++] = RTYSC14A_IKrInit; 
      }
      cInit[k++] = OHaraRudy_FluxesInit;      // must be initialize after ICa
      assert(k == nComp); 
      privateStateOffset_=(int *)malloc(sizeof(int)*nComp); 
      privateParmsOffset_=(int *)malloc(sizeof(int)*nComp); 
      info_=(COMPONENTINFO *)malloc(sizeof(COMPONENTINFO)*nComp); 
      int stateOffset =0; 
      int parmsOffset =0; 
      stateSize_=0; 
      parmsSize_=0; 
      k=0; 
      for (int i=0;i<nComp;i++)
      {
         privateStateOffset_[i] =  stateOffset;
         privateParmsOffset_[i] =  parmsOffset;
         COMPONENTINFO info = cInit[i](); 
         info_[i] = info; 
         VARINFO *varInfo = info.varInfo; 
         for (int j=0;j<info.nVar;j++) 
         {
            string name = varInfo[j].name ;
            string units = varInfo[j].units;
            unsigned int index = (i << 16) + varInfo[j].index;
            bool  checkpoint = false;
            if (varInfo[j].type == PSTATE_TYPE)
            {
               stateOffset++;
               if (  name != "Vm" ) checkpoint = true; 
               else 
               {
                  if (cellType == ENDO_CELL) defaultVoltage_ = varInfo[j].defaultValueENDO; 
                  if (cellType == EPI_CELL) defaultVoltage_ = varInfo[j].defaultValueEPI; 
                  if (cellType == M_CELL) defaultVoltage_ = varInfo[j].defaultValueM; 
               }
            }
            if (varInfo[j].type == PARAMETER_TYPE)
            {
               parmsOffset++;
            }
            handleMap_[name] = CheckpointVarInfo(index, checkpoint, units );
         }
      } 
      printf("size = %d %d\n",stateOffset,parmsOffset); 
      stateSize_ = stateOffset*sizeof(double); 
      parmsSize_ = parmsOffset*sizeof(double); 

      char **nameState = (char **)malloc(stateSize_/8*sizeof(char *)); 
      defaultStateENDO_ = (double *)malloc(stateSize_); 
      defaultStateEPI_ = (double *)malloc(stateSize_); 
      defaultStateM_ = (double *)malloc(stateSize_); 

      char **nameParms = (char **)malloc(parmsSize_/8*sizeof(char *)); 
      defaultParmsENDO_ = (double *)malloc(parmsSize_); 
      defaultParmsEPI_  = (double *)malloc(parmsSize_); 
      defaultParmsM_    = (double *)malloc(parmsSize_); 
      int nState = 0; 
      int nParms = 0; 
      for (int i=0;i<nComp;i++)
      {
         COMPONENTINFO info = info_[i]; 
         VARINFO *varInfo = info.varInfo; 
         for (int j=0;j<info.nVar;j++) 
         {
            string name = varInfo[j].name ;
            double defaultValueENDO  = varInfo[j].defaultValueENDO;
            double defaultValueEPI = varInfo[j].defaultValueEPI  ;
            double defaultValueM  = varInfo[j].defaultValueM;
            if (varInfo[j].type == PSTATE_TYPE)
            {
               nameState[nState] = strdup(name.c_str()); 
               defaultStateENDO_[nState] = defaultValueENDO; 
               defaultStateEPI_[nState] = defaultValueEPI; 
               defaultStateM_[nState] = defaultValueM; 
               nState++; 
            }
            if (varInfo[j].type == PARAMETER_TYPE)
            {
               nameParms[nParms] = strdup(name.c_str()); 
               defaultParmsENDO_[nParms] = defaultValueENDO; 
               defaultParmsEPI_[nParms] = defaultValueEPI; 
               defaultParmsM_[nParms] = defaultValueM; 
               nParms++; 
            }
         }
      }
      assert(nState*sizeof(double) == stateSize_); 
      initialized = true; 
   for (int i=0;i<nState;i++) printf("OHaraRudy %2d %16s %16.8e %16.8e %16.8e\n",i,nameState[i],defaultStateENDO_[i],defaultStateEPI_[i],defaultStateM_[i]); 
   for (int i=0;i<nParms;i++) printf("OHaraRudy %2d %16s %16.8e %16.8e %16.8e\n",i,nameParms[i],defaultParmsENDO_[i],defaultParmsEPI_[i],defaultParmsM_[i]); 
   }
}

/* Remember that down in the cell models the units don't necessarily
 *  correspond to the internal units of Cardioid.  The units in this map
 *  are the units the cell model expects the variables to have.
 */
HandleMap& OHaraRudy::getHandleMap()
{
   return handleMap_;
}

void OHaraRudy::initState(int cellType)
{
   int nState = stateSize_/sizeof(double); 
   if (cellType == ENDO_CELL) for (int i=0;i<nState;i++) state_[i] = defaultStateENDO_[i]; 
   if (cellType == EPI_CELL) for (int i=0;i<nState;i++)  state_[i] = defaultStateEPI_[i]; 
   if (cellType == M_CELL) for (int i=0;i<nState;i++)    state_[i] = defaultStateM_[i]; 
}
void OHaraRudy::initParms(int cellType)
{
   if (cellType == ENDO_CELL)  cellParms_ = defaultParmsENDO_; 
   if (cellType ==  EPI_CELL)  cellParms_ = defaultParmsEPI_; 
   if (cellType ==    M_CELL)  cellParms_ = defaultParmsM_;  
}
void OHaraRudy::setValue(int varHandle, double value) 
{
   int iComp = varHandle >> 16; 
   int index = varHandle &  0xffff; 
   double *pstate = state_ + privateStateOffset_[iComp]; 
   double *cellParms = cellParms_ + privateParmsOffset_[iComp]; 
   info_[iComp].access(WRITE,index,&value,cellParms,pstate); 
}
double OHaraRudy::getValue(int varHandle)  const
{
   double value; 
   int iComp = varHandle >> 16; 
   int index = varHandle  & 0xffff; 
   double *pstate = state_ + privateStateOffset_[iComp]; 
   double *cellParms = cellParms_ + privateParmsOffset_[iComp]; 
   info_[iComp].access(READ,index,&value,cellParms,pstate); 
   return value; 
}
