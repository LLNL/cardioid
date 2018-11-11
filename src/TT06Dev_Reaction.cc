#include "TT06Dev_Reaction.hh"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "Anatomy.hh"
#include "TT06Func.hh"
#include "TT06Gates.h"
#include "TT06NonGates.h"
#include "pade.hh"
#include "PerformanceTimers.hh"
#include "ThreadServer.hh"
#include "mpiUtils.h"
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <algorithm>
#include "pio.h"
#include "object_cc.hh"
#include "ioUtils.h"
#include "reactionFactory.hh"

#include "slow_fix.hh"

//#define SHOWVAR(x) std::cout << #x " = " << x << std::endl


using namespace std;
using namespace PerformanceTimers;
using namespace TT06Func;
int workBundle(int index, int nItems, int nGroups , int mod, int& offset);

#if 0
static const  UPDATEGATE gateEq[] ={ update_mGate_v1, update_hGate_v1, update_jGate_v1, update_Xr1Gate_v1, 
                               update_Xr2Gate_v1, update_XsGate_v1, update_rGate_v1, update_dGate_v1, 
                               update_fGate_v1, update_f2Gate_v1,  update_jLGate_v1, update_s0Gate_v1, 
                               update_s1Gate_v1} ;
#else
static const  UPDATEGATE gateEq[] ={ update_mGate_v2, update_hGate_v2, update_jGate_v2, update_Xr1Gate_v2, 
                               update_Xr2Gate_v2, update_XsGate_v2, update_rGate_v2, update_dGate_v2, 
                               update_fGate_v2, update_f2Gate_v2,  update_jLGate_v2, update_s0Gate_v2, 
                               update_s1Gate_v2} ;
#endif


OVF fitFuncMap(string name) 
{
      const char *fitName[] ={ "fv0", "fv1", "fv2", "fv3", "fv4", "fv5", "fv6", "mMhu", "mTauR", "hMhu", "hTauR", "hTauRMod", "jMhu", "jTauR", "jTauRMod", "Xr1Mhu", "Xr1TauR", "Xr2Mhu", "Xr2TauR", 
                           "XsMhu", "XsTauR", "rMhu", "rTauR", "dMhu", "dTauR", "fMhu", "fTauR", "f2Mhu", "f2TauR", "jLMhu", "jLTauR", "sMhu0", "sTauR0", "sMhu1", "sTauR1" ,NULL}; 
      OVF func[] =      { fv0 ,  fv1 ,  fv2 ,  fv3 ,  fv4 ,  fv5 , fv6 ,  mMhu ,   mTauR ,  hjMhu , hTauR ,  hTauRMod , hjMhu ,  jTauR ,  jTauRMod ,  Xr1Mhu ,  Xr1TauR ,  Xr2Mhu ,  Xr2TauR , 
                            XsMhu ,  XsTauR ,  rMhu ,  rTauR ,  dMhu ,  dTauR ,  fMhu ,  fTauR ,  f2Mhu ,  f2TauR ,  jLMhu ,  jLTauR ,  sMhu0 ,  sTauR0 ,  sMhu1 ,  sTauR1 , NULL}; 
      int i=0; 
      while ( fitName[i] != NULL) 
      {
	if (strcmp(fitName[i],name.c_str())==0)  return func[i]; 
        i++; 
      }
      return NULL; 
}
void writeFit(string fitFile, PADE *fit, int cnt)
{
   if (getRank(0)==0) 
   {
      FILE *file=fopen(fitFile.c_str(),"w"); 
      fprintf(file,"functions FIT { functions="); 
      for (int index =0;index< cnt;index++)
           fprintf(file," %s" ,  fit[index].name.c_str()); 
      fprintf(file,";}\n"); 
      for (int index =0;index< cnt;index++)
      
      {
         padeErrorInfo(fit[index],index); 
         padeWrite(file,fit[index]); 
      }
      fclose(file); 
   }
}

// Can't pass parms by const reference since the parms contain a map and
// we access the map with the subscript operator.  This is a non-const
// operation.  Maybe we can re-write this code to avoid the subscript
// operator. 
TT06Dev_Reaction::TT06Dev_Reaction(const double dt, const int numPoints, TT06Dev_ReactionParms& parms, const ThreadTeam& group)
                 : nCells_(numPoints), group_(group)
{
   static bool initialized = false;
   if (! initialized)
   {
      initialized = true;
      TT06Func::initCnst();
      initExp(); 
      int pid;
      MPI_Comm_rank(MPI_COMM_WORLD,&pid);	  
      string fCassFormName[] = { "TT06", "RICE"}; 
      if (pid ==0) printf("fCassForm = %s\n",fCassFormName[fCassForm].c_str()); 
   }   
   nCellBuffer_ =  convertActualSizeToBufferSize(nCells_);
   mkCellTypeParms_(parms);
   for (int ii=0; ii<11; ii++) {
      gateThreadMap_[ii] = ii;
   }
   if (cellTypeParm_.s_switch) {
      gateThreadMap_[11] = 12;
   } else {
      gateThreadMap_[11] = 11;
   }

   mkState_(parms);

   fastReaction_ = 0; 
   //if (parms.fastReaction >= 0) ;
   {
      if ( (fabs(parms.tolerance/1e-4 - 1.0) < 1e-12) && parms.jhTauSmooth == 1 )fastReaction_ = parms.fastReaction; 
   }
   int fastGates =          fastReaction_&0x00ff; 
   if ( fastGates   == 0) update_gate_   =TT06Func::updateGate;
   if ( fastGates   >= 1) update_gate_   =TT06Func::updateGateFast;

   int fastNonGates = fastReaction_/256 & 0x00ff; 
   if (fastNonGates  == 0) update_nonGate_=update_nonGate;
   if (fastNonGates  == 1) update_nonGate_=update_nonGateSimdM;
   if (fastNonGates  == 2) update_nonGate_=update_nonGateSimdF;
   if (fastNonGates  >= 3) update_nonGate_=update_nonGateSimdFA;

   int ompId =  omp_get_thread_num(); 
   //const int tid = group_.teamRank();

   mkFitParms_(parms);

   PADE *gateFit = fit_+gateFitOffset; 
   for (int eq=0;eq<13;eq++) 
   {
      double *tauR  = gateFit[2*eq+1].coef; 
      int m         = gateFit[2*eq+1].m; 
      for (int j=0;j<m;j++)    tauR[j] *= dt;
   }

   double **gate = &state_[gateOffset]; 
   int eq; 
   for (eq=0;eq<11;eq++)
   {
      gateEqX_[eq]  = gateEq[eq];  
      gateX_[eq]    = gate[eq];
      mhuX_[eq]     = gateFit[2*eq+0].coef; 
      tauRX_[eq]    = gateFit[2*eq+1].coef; 
   }
   //sGate
   for (eq=11;eq<13;eq++)
   {
      gateEqX_[eq]  = gateEq[eq];
      gateX_[eq] = gate[11];
      mhuX_[eq]  = gateFit[2*eq+0].coef;
      tauRX_[eq] = gateFit[2*eq+1].coef;
   }

   mkWorkBundles_(parms);
}
void TT06Dev_Reaction::mkCellTypeParms_(TT06Dev_ReactionParms& parms)
{
   cellTypeParm_.s_switch = parms.cellTypeParm.s_switch;  
   cellTypeParm_.P_NaK    = parms.cellTypeParm.P_NaK;  
   cellTypeParm_.g_Ks     = parms.cellTypeParm.g_Ks ;  
   cellTypeParm_.g_Kr     = parms.cellTypeParm.g_Kr ;  
   cellTypeParm_.g_to     = parms.cellTypeParm.g_to ;  
   cellTypeParm_.g_NaL    = parms.cellTypeParm.g_NaL;  
   XXXinitialVm_ = parms.cellTypeParm.Vm;
   XXXstateInitial_.resize(nStateVar);
   int stateCnt=0; 
   map<string,STATE>stateMap = parms.cellTypeParm.state; 
   for (map<string,STATE>::iterator it=stateMap.begin();it!=stateMap.end();it++)
   {
      STATE stateDefault = it->second; 
      double value = stateDefault.value; 
      int type = stateDefault.type; 
      int index = stateDefault.index; 
      XXXstateInitial_[index]=value; stateCnt++;
   }
   assert(stateCnt==nStateVar); 
}
void TT06Dev_Reaction::mkState_(TT06Dev_ReactionParms& parms)
{
   unsigned bufSize = sizeof(double)*nStateVar*nCellBuffer_;
   int rc = posix_memalign((void**)&stateBuffer_, 32, bufSize);
   assert((size_t)stateBuffer_ % 32 == 0);


   state_.resize(nStateVar); 
   for (unsigned jj=0; jj<nStateVar; ++jj)
   {
      state_[jj]= stateBuffer_+jj*nCellBuffer_;
      assert((size_t)state_[jj] % 32 == 0);
   }

   double c9=get_c9(); 
   for (unsigned ii=0; ii<nCellBuffer_; ++ii)
   {
      for (int j=0;j<nStateVar;j++) state_[j][ii]  = XXXstateInitial_[j]; 
      state_[dVK_i][ii] = state_[K_i][ii]/c9+XXXinitialVm_;
   }
}
void TT06Dev_Reaction::mkFitParms_(TT06Dev_ReactionParms& parms)
{
   fit_ = parms.fit; 
   if (fit_ == NULL ) 
   {
      const char *fitName[] ={ "fv0", "fv1", "fv2", "fv3", "fv4", "fv5", "fv6", "mMhu", "mTauR", "hMhu", "hTauR",  "jMhu", "jTauR", "Xr1Mhu", "Xr1TauR", "Xr2Mhu", "Xr2TauR", 
         "XsMhu", "XsTauR", "rMhu", "rTauR", "dMhu", "dTauR", "fMhu", "fTauR", "f2Mhu", "f2TauR", "jLMhu", "jLTauR", "sMhu0", "sTauR0", "sMhu1", "sTauR1" ,NULL}; 
      int cnt=0; 
      for ( int i=0; fitName[i] != NULL;i++) cnt++; 
      fit_ =  new PADE  [cnt] ;
      int maxCost=128; 
      int maxTerms=64; 
      double tol = parms.tolerance; 
      double deltaV = 0.1; 
      for ( int i=0; i < cnt; i++) 
      {
         string name = fitName[i]; 
         if (parms.jhTauSmooth  && i == 10) name = "hTauRMod";
         if (parms.jhTauSmooth  && i == 12) name = "jTauRMod";
         double V0 = (i == 6) ?   0 : -100.0; 
         double V1 = (i == 6) ? 130 :   50.0; 
         padeApprox(fit_[i],name,fitFuncMap(name),NULL,0,deltaV,V0,V1,tol,maxTerms,maxTerms,maxCost,0,0,NULL); 
         //padeSet(fit_+i,maxTerms,maxTerms,maxCost); 
      }
      writeFit(parms.fitFile,fit_,cnt); 
   }
}
void TT06Dev_Reaction::mkWorkBundles_(TT06Dev_ReactionParms& parms)
{

   int nThreads = group_.nThreads();
   int nSquads =  group_.nSquads();
 if (nThreads ==0 && nSquads ==0) return;
   assert(nThreads%nSquads==0) ;
   int squadSize=1; 
   if (nSquads >0) squadSize= nThreads/nSquads;
   assert(12 % squadSize == 0);
   int nEq = 12/squadSize;
   int sGateIndex = 11; 
   if (parms.gateThreadMap.size() > 0) 
   {
      assert(parms.gateThreadMap.size()  == 12); 
      for (int ii=0;ii<12;ii++) 
      {
         int m = parms.gateThreadMap[ii];
         if ( m == 11 && cellTypeParm_.s_switch) {
            m = 12;
         }
         gateThreadMap_[ii] = m;
      }
   }
   gateWork_.resize(nThreads); 
   nonGateWork_.resize(nThreads);   
   int splitRank =-1; 
   int nOmpThreads =  omp_get_max_threads();
   for (int id=0;id<nOmpThreads;id++) 
   {

      const ThreadRankInfo& rankInfo = group_.rankInfo(id);
      int teamRank = rankInfo.teamRank_; 
      if (teamRank == -1) continue; 
      int squadId   =rankInfo.coreRank_; 
      int squadRank =rankInfo.squadRank_; 
      assert(squadSize == rankInfo.squadSize_);
      int offset; 
      int nCell = workBundle(squadId, nCells_, nSquads , 4, offset);
      int offsetEq = nEq*squadRank; 
      gateWork_[teamRank].offsetCell =  offset; 
      gateWork_[teamRank].nCell    =  nCell; 
      gateWork_[teamRank].nEq      = nEq; 
      gateWork_[teamRank].map      = gateThreadMap_+offsetEq; 

      /* non Gate work partitioning */ 
      int nv = (nCell+3)/4;
      int q = nv/squadSize; 
      int r = nv%squadSize;
      int i0,i1;

      i0 = q*squadRank;
      if(squadRank < r) 
      {
         i0 = i0 + squadRank;
         i1 = i0 + q + 1;
      }       
      else 
      {
         i0 = i0 + r;
         i1 = i0 + q;
      }
      i0 = i0*4;
      i1 = i1*4;
      if(i1 > nCell) i1 = nCell;
      nonGateWork_[teamRank].offsetCell = offset + i0;
      nonGateWork_[teamRank].nCell = i1 - i0;
   }
}
TT06Dev_Reaction::~TT06Dev_Reaction()
{
   free(stateBuffer_);
   delete fit_;
}
// void TT06Dev_Reaction::writeStateDev(int loop)
// {
//    int  map[] = { dVK_i , Na_i , Ca_i , Xr1_gate , Xr2_gate , Xs_gate , m_gate , h_gate , j_gate , Ca_ss , d_gate , f_gate , f2_gate , fCass , s_gate , r_gate , Ca_SR , R_prime , jL_gate};
//    int  isGate[] = {    0  ,     0, 0 , 1 , 1 , 1 , 1 , 1 , 1 , 0 , 1 , 1 , 1 , 0 , 1 , 1 , 0 , 0 , 1};
//    for (int i=0;i<nStateVar;i++) 
//    {
//       int k = map[i]; 
//       double state; 
//       state = state_[k][0];
//       printf("%d %24.14le\n",i,state); 
//    }
// }


int workBundle(int index, int nItems, int nGroups, int mod, int& offset)
{
   assert(0<=index && index < nGroups); 
   int nItems0 = (nItems+mod-1)/mod ;
   int n = nItems0/nGroups;
   int remainder = nItems0-n*nGroups;
   int m = nGroups - remainder;
   offset =0;
   if ( index <  m )
   {
      offset += n*index;
   }
   else
   {
      offset += n*m + (n+1)*(index-m);
      n++;
   }
   n      *= mod;
   offset *= mod;
   //if ( index == nGroups -1) n = mod*((nItems-offset)/mod);
   if ( index == nGroups -1) n = ((nItems-offset));

   assert (offset%4 == 0);

   return n;
}

int partition(int index, int nItems, int nGroups , int& offset)
{
   int n = nItems/nGroups; 
   int remainder = nItems-n*nGroups; 
   int m = nGroups - remainder; 
   if ( index <  m ) 
   {
      offset += n*index; 
   }
   else
   {
      offset += n*m + (n+1)*(index-m); 
      n++; 
   }
   return n; 
}
int TT06Dev_Reaction::nonGateWorkPartition_(int& offset)
{
   offset=0; 

   const int simdLength = 4;
   int nSimd = nCells_/simdLength;
   int simdRemainder = nCells_%simdLength;

   int nThreads = group_.nThreads();
   int nSimdPerThread = nSimd / nThreads;
   int threadRemainder = nSimd % nThreads;

   int maxAssignableThreads = nSimd + (simdRemainder > 0);
   int nActiveThreads = min(nThreads, maxAssignableThreads);

   int threadRank = group_.teamRank();
   if (threadRank >= nActiveThreads)
   {
      // too few cells to keep all threads busy.
      offset = 0;
      return 0;
   }
   if (threadRank < threadRemainder)
   {
      offset = (nSimdPerThread + 1) * threadRank * simdLength;
      return (nSimdPerThread + 1) * simdLength;
   }
   offset = (threadRemainder + (nSimdPerThread * threadRank)) * simdLength;
   int size = nSimdPerThread * simdLength;
   if (threadRank + 1 == nActiveThreads)
      size += simdRemainder;
   return size;
}

void TT06Dev_Reaction::calc(double dt,
                            ro_mgarray_ptr<double> Vm_m,
                            ro_mgarray_ptr<double> iStim_m,
                            wo_mgarray_ptr<double> dVm_m)
{
   ro_array_ptr<double> Vm = Vm_m.useOn(CPU);
   ro_array_ptr<double> iStim = iStim_m.useOn(CPU);
   wo_array_ptr<double> dVm = dVm_m.useOn(CPU);
   WORK work ={ 0,static_cast<int>(nCells_),0,12}; 
   if (nCells_ > 0) 
   {
      update_nonGate_((void*)fit_,dt,&(cellTypeParm_), nCells_,const_cast<double *>(&Vm[0]),  0, &state_[0], dVm.raw());
      update_gate_(dt, nCells_, cellTypeParm_.s_switch, const_cast<double *>(&Vm[0]), 0, &state_[gateOffset],fit_,work);
   }
}
void TT06Dev_Reaction::updateNonGate(double dt, ro_mgarray_ptr<double> Vm_m, wo_mgarray_ptr<double> dVR_m)
{
   int offset,nCells; 

   ro_array_ptr<double> Vm = Vm_m.useOn(CPU);
   wo_array_ptr<double> dVR = dVR_m.useOn(CPU);
   
#ifdef LEGACY_NG_WORKPARTITION
   nCells = nonGateWorkPartition(offset);
#else
   int ompId =  omp_get_thread_num(); 
   const int tid = group_.teamRank();
   offset = nonGateWork_[tid].offsetCell;
   nCells = nonGateWork_[tid].nCell;
#endif

   if (nCells > 0) 
   {
      startTimer(nonGateTimer);
      update_nonGate_(fit_,dt,&cellTypeParm_, nCells, const_cast<double *>(&Vm[offset]),  offset, &state_[0], dVR.raw()+offset);
      stopTimer(nonGateTimer);
   }
}
void TT06Dev_Reaction::updateGate   (double dt, ro_mgarray_ptr<double> Vm_m)
{
   int teamRank = group_.teamRank();

   ro_array_ptr<double> Vm = Vm_m.useOn(CPU);

   TT06Func::WORK work = gateWork_[teamRank]; 
   int nCell=work.nCell; 
   if ( nCell ==0)  return; 
   int offsetCell=work.offsetCell; 
   int *map = work.map; 
   int nEq  = work.nEq;       //We could assume each HW thread get the same number of equations; 

   double *vm = const_cast<double *>(&Vm[offsetCell]); 

   startTimer(gateTimer);
   for (int ii=0;ii<nEq;ii++)
   {
      int eq = map[ii]; 
      gateEqX_[eq](dt, nCell , vm , gateX_[eq]+offsetCell, mhuX_[eq], tauRX_[eq]);
   }
   stopTimer(gateTimer);
}

void TT06Dev_Reaction::initializeMembraneVoltage(wo_mgarray_ptr<double> Vm_m)
{
   assert(Vm_m.size() >= nCells_);
   wo_array_ptr<double> Vm = Vm_m.useOn(CPU);
   for (unsigned ii=0; ii<nCells_; ++ii) Vm[ii] = XXXinitialVm_;
}

void TT06Dev_Reaction::getCheckpointInfo(vector<string>& name,
      vector<string>& unit) const
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

/** This function maps the string representation of a variable name to
 *  the handle representation.  Returns the value -1 for
 *  unrecognized varName. */
int TT06Dev_Reaction::getVarHandle(const string& varName) const
{
   return getHandleMap()[varName].handle_;
}

void TT06Dev_Reaction::setValue(int iCell, int varHandle, double value)
{
   assert(varHandle >= 0);

   if (varHandle < nStateVar)
      state_[varHandle][iCell] = value;

   // no support for variables other than state variables at present.
   assert(varHandle < nStateVar);

}

double TT06Dev_Reaction::getValue(int iCell, int handle) const
{
   assert(handle >= 0);

   if (handle < nStateVar)
      return state_[handle][iCell];

   // no support for variables other than state variables at present.
   assert(handle < nStateVar);

   return 0.;
}

void TT06Dev_Reaction::getValue(int iCell,
      const vector<int>& handle,
      vector<double>& value) const
{
   for (unsigned ii=0; ii<handle.size(); ++ii)
      value[ii] = getValue(iCell, handle[ii]);
}


const string TT06Dev_Reaction::getUnit(const string& varName) const
{
   return getHandleMap()[varName].unit_;
}


/** Remember that down in the cell models the units don't necessarily
 *  correspond to the internal units of Cardioid.  The units in this map
 *  are the units the cell model expects the variables to have. */
HandleMap& TT06Dev_Reaction::getHandleMap() const
{
   static HandleMap handleMap;
   if (handleMap.size() == 0)
   {
      handleMap["dVK_i"]      = CheckpointVarInfo(dVK_i,      true,  "mV");
      handleMap["Na_i"]       = CheckpointVarInfo(Na_i,       true,  "mM");
      handleMap["Ca_i"]       = CheckpointVarInfo(Ca_i,       true,  "mM");
      handleMap["Xr1_gate"]   = CheckpointVarInfo(Xr1_gate,   true,  "1");
      handleMap["Xr2_gate"]   = CheckpointVarInfo(Xr2_gate,   true,  "1");
      handleMap["Xs_gate"]    = CheckpointVarInfo(Xs_gate,    true,  "1");
      handleMap["m_gate"]     = CheckpointVarInfo(m_gate,     true,  "1");
      handleMap["h_gate"]     = CheckpointVarInfo(h_gate,     true,  "1");
      handleMap["j_gate"]     = CheckpointVarInfo(j_gate,     true,  "1");
      handleMap["Ca_ss"]      = CheckpointVarInfo(Ca_ss,      true,  "mM");
      handleMap["d_gate"]     = CheckpointVarInfo(d_gate,     true,  "1");
      handleMap["f_gate"]     = CheckpointVarInfo(f_gate,     true,  "1");
      handleMap["f2_gate"]    = CheckpointVarInfo(f2_gate,    true,  "1");
      handleMap["fCass_gate"] = CheckpointVarInfo(fCass,      true,  "1");
      handleMap["s_gate"]     = CheckpointVarInfo(s_gate,     true,  "1");
      handleMap["r_gate"]     = CheckpointVarInfo(r_gate,     true,  "1");
      handleMap["Ca_SR"]      = CheckpointVarInfo(Ca_SR,      true,  "mM");
      handleMap["R_prime"]    = CheckpointVarInfo(R_prime,    true,  "1");
      handleMap["jL_gate"]    = CheckpointVarInfo(jL_gate,    true,  "mM");
      assert(handleMap.size() == nStateVar);
   }
   return handleMap;
}

   PADE *scanFit()
   {
      if (!object_exists("functions", "FIT") ) return NULL;
      OBJECT* obj = object_find("functions", "FIT");
      vector<string>fitName; 
      objectGet(obj,"functions",fitName); 
      PADE *fit =  new PADE  [fitName.size()] ;
      for (int i=0;i<fitName.size();i++)
      {
         vector<double> coef; 
         double tol,deltaV,V0,V1;
         int l,m;
         string name = fitName[i]; 
         OBJECT* obj = object_find(name.c_str(), "FIT");
         objectGet(obj,"tol", tol,"1.0"); 
         objectGet(obj,"deltaX",deltaV,"1"); 
         objectGet(obj,"x0"    ,    V0,"0"); 
         objectGet(obj,"x1"    ,    V1,"0"); 
         objectGet(obj,"l"     ,     l,"0"); 
         objectGet(obj,"m"     ,     m,"0"); 
         objectGet(obj,"coef"  ,coef); 
         padeApprox(fit[i],name,fitFuncMap(name),NULL,0,deltaV,V0,V1,tol,0,0,0,l,m,&coef[0]); 
      //   fit[i].aparms=fit+i; 
      }
      return fit; 
   }
   REACTION_FACTORY(TT06Dev)(OBJECT* obj, const double dt, const int numPoints, const ThreadTeam& group)
   {
      vector<string> scaleCurrents;
      TT06Dev_ReactionParms parms;
      map<string,TT06Func::CellTypeParmsFull> cellTypeParms = TT06Func::getStandardCellTypes(); 
      int fastGate =-1; 
      int fastNonGate =-1; 
      TT06Func::initCnst(); 
      string fitFit; 
/*
      if (object_testforkeyword(obj,"mod") )
      {
          //if (getRank(0) == 0) printf("keyword <mod> is deprecated by <jhTauSmooth>.  Replace keyword <mod> in the reaction object <%s>  with <jhTauSmooth>\n",obj->name); 
          
          //exit(1); 
         objectGet(obj, "mod",          parms.jhTauSmooth, "0") ;
      }
*/
      objectGet(obj, "mod",          parms.jhTauSmooth, "0") ;
      objectGet(obj, "tolerance",    parms.tolerance, "0.0") ;
      objectGet(obj, "fastReaction", parms.fastReaction, "-1") ;
      objectGet(obj, "fastGate",     fastGate, "-1") ;
      objectGet(obj, "fastNonGate",  fastNonGate, "-1") ;
      if (object_testforkeyword(obj, "gateThreadMap") ) objectGet(obj,"gateThreadMap",parms.gateThreadMap);
      parms.fit = NULL; 
      parms.fitFile = "fit.data"; 
      if (object_testforkeyword(obj,"fitFile") )
      {
        objectGet(obj, "fitFile",     parms.fitFile,"fit.data"); 
        int fileFitExists=0; 
        if (getRank(0) ==0) 
        {
          if (filetest(parms.fitFile.c_str(), S_IFREG) == 0) fileFitExists=1; 
          if (fileFitExists) object_compilefile(parms.fitFile.c_str()); 
        }
        MPI_Bcast(&fileFitExists, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (fileFitExists) 
        {
           object_Bcast(0, MPI_COMM_WORLD);
           parms.fit = scanFit();
        }
      }
      if (parms.fastReaction == -1) 
      {
         if (fastGate    ==  -1 )  fastGate   =0; 
         if (fastNonGate ==  -1 )  fastNonGate=0; 
         parms.fastReaction = fastGate+256*fastNonGate; 
      }
      string cellTypeName;
      objectGet(obj, "cellTypeName", cellTypeName, "endoCellML");

      OBJECT* cellobj = object_find2(cellTypeName.c_str(), "CELLTYPE", IGNORE_IF_NOT_FOUND);
      if (cellobj) {
         string clone; 
         objectGet(cellobj, "clone", clone, "") ;
         if (clone != "")
         {
            assert(cellTypeName != clone); 
            assert(cellTypeParms.count(clone) != 0);
            cellTypeParms[cellTypeName]=cellTypeParms[clone]; 
            cellTypeParms[cellTypeName].name=cellTypeName;
         }
         if (object_testforkeyword(cellobj, "s_switch")       ) objectGet(cellobj,"s_switch"      ,cellTypeParms[cellTypeName].s_switch,"0"); 
         if (object_testforkeyword(cellobj, "P_NaK")          ) objectGet(cellobj,"P_NaK"         ,cellTypeParms[cellTypeName].P_NaK,"0.0"); 
         if (object_testforkeyword(cellobj, "g_NaL")          ) objectGet(cellobj,"g_NaL"         ,cellTypeParms[cellTypeName].g_NaL,"0.0"); 
         if (object_testforkeyword(cellobj, "g_Ks")           ) objectGet(cellobj,"g_Ks"        ,cellTypeParms[cellTypeName].g_Ks,"0.0"); 
         if (object_testforkeyword(cellobj, "g_Kr")           ) objectGet(cellobj,"g_Kr"        ,cellTypeParms[cellTypeName].g_Kr,"0.0"); 
         if (object_testforkeyword(cellobj, "g_to")           ) objectGet(cellobj,"g_to"        ,cellTypeParms[cellTypeName].g_to,"0.0"); 
      } else {
         assert(cellTypeName == cellTypeParms[cellTypeName].name);
      }
      parms.cellTypeParm = cellTypeParms[cellTypeName];

      parms.currentNames = scaleCurrents;
      
      Reaction *reaction = new TT06Dev_Reaction(dt, numPoints, parms, group);
      return  reaction; 
   }
   REACTION_FACTORY(TT06Opt)(OBJECT* obj, const double dt, const int numPoints, const ThreadTeam& group)
   {
      return reactionFactoryForTT06Dev(obj,dt,numPoints,group);
   }
