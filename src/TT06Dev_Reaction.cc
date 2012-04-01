#include "TT06Dev_Reaction.hh"
#include <cmath>
#include "Anatomy.hh"
#include "TT06Func.hh"
#include "pade.hh"
#include "PerformanceTimers.hh"
#include "ThreadServer.hh"
#include <cstdlib>

using namespace std;
using namespace TT06Func;
using namespace PerformanceTimers;

static FILE *tfile[64]; 

TT06Dev_Reaction::TT06Dev_Reaction(const Anatomy& anatomy,double tolerance,int mod, const ThreadTeam& group)
: nCells_(anatomy.nLocal()),
  group_(group)
{
   static bool initialized = false;
   if (! initialized)
   {
      initialized = true;
      TT06Func::initCnst();
   }
   {
      double V0 = -100.0; 
      double V1 =  50.0; 
      double deltaV = 0.1; 
      int maxCost=128; 
      int maxTerms=64; 
      PADE **fit=TT06Func::makeFit(tolerance,V0,V1,deltaV,mod); 
      for (int i=0;fit[i]!=NULL;i++) padeCalc(fit[i],maxTerms,maxTerms,maxCost); 
      TT06Func::writeFit(fit); 
   }
   dtForFit_=0.0; 
   ttType_.resize(256, -1); 
   ttType_[30] = 0;
   ttType_[31] = 0;
   ttType_[75] = 0;
   ttType_[76] = 1;
   ttType_[77] = 2;
   ttType_[100] = 0;
   ttType_[101] = 1;
   ttType_[102] = 2;
   s_.resize(nCells_);
   gates_ = (double **) malloc(sizeof(double *)*nGatesVar); 
   for (unsigned jj=0; jj<nGatesVar; ++jj) gates_[jj]= (double *)malloc(sizeof(double)*nCells_); 
   for (unsigned ii=0; ii<nCells_; ++ii)
   {
      double gate[nGatesVar]; 
      assert(anatomy.cellType(ii) >= 0 && anatomy.cellType(ii) < 256);
      int cellType = ttType_[anatomy.cellType(ii)];
      TT06Func::initState(&(s_[ii]),gate,cellType); 
      TT06Func::initGate(gates_,ii,gate); 
   }
   
}
TT06Dev_Reaction::~TT06Dev_Reaction()
{
}
void TT06Dev_Reaction::writeStateDev(int loop)
{
   int  map[] = { dVK_i , Na_i , Ca_i , Xr1_gateN , Xr2_gateN , Xs_gateN , m_gateN , h_gateN , j_gateN , Ca_ss , d_gateN , f_gateN , f2_gateN , fCass_gate , s_gateN , r_gateN , Ca_SR , R_prime , jL_gateN};
   int  isGate[] = {    0  ,     0, 0 , 1 , 1 , 1 , 1 , 1 , 1 , 0 , 1 , 1 , 1 , 0 , 1 , 1 , 0 , 0 , 1};
   for (int i=0;i<nStateVar+nGatesVar;i++) 
   {
      int k = map[i]; 
      double state; 
      if (isGate[i]) 
      {
         state = gates_[k][0];
      }
      else 
      {
         state = s_[0].state[k];
      }
      printf("%d %24.14le\n",i,state); 
   }
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
int TT06Dev_Reaction::nonGateWorkPartition(int& offset)
{
//     int coreID,hwThreadID,threadID,nCores,nHwThreads,nThreads;
//     group_.groupInfo(coreID,hwThreadID,threadID,nCores,nHwThreads,nThreads); 
   offset=0; 
   
   const ThreadRankInfo& rankInfo = group_.rankInfo();
   int nCores = group_.nSquads();
   int coreID = rankInfo.coreRank_;
   int hwThreadID = rankInfo.squadRank_;
   int nHwThreads = rankInfo.squadSize_;
   int nCellsCore= partition(coreID, nCells_, nCores, offset);
   int nCells = partition(hwThreadID, nCellsCore, nHwThreads, offset); 
   return nCells; 
}

void TT06Dev_Reaction::calc(double dt, const vector<double>& Vm, const vector<double>& iStim, vector<double>& dVm)
{
   TT06Func::updateNonGate(dt, nCells_,&Vm[0], &(s_[0]), 0, gates_, &dVm[0]);
   TT06Func::updateGate(dt, nCells_,&Vm[0], &(s_[0]),0,gates_);
}
void TT06Dev_Reaction::updateNonGate(double dt, const vector<double>& Vm, vector<double>& dVR)
{
   int offset; 
   profileStart(nonGateTimer);
   int nCells = nonGateWorkPartition(offset); 
   TT06Func::updateNonGate(dt, nCells,&Vm[offset], &(s_[offset]), offset, gates_, &dVR[offset]);
   profileStop(nonGateTimer);
}
void TT06Dev_Reaction::updateGate(double dt, const vector<double>& Vm)
{
   int offset=0; 
   int nCells=nCells_; 
   profileStart(gateTimer);
#if (1) 
   nCells = nonGateWorkPartition(offset); 
   TT06Func::updateGate(dt, nCells,&Vm[offset], &(s_[offset]),offset,gates_);
#else
   int id = groupThreadID(group_); 
   if (id ==0) TT06Func::updateGate0(dt, nCells,&Vm[offset], &(s_[offset]),offset,gates_);
   if (id ==1) TT06Func::updateGate1(dt, nCells,&Vm[offset], &(s_[offset]),offset,gates_);
   if (id ==2) TT06Func::updateGate2(dt, nCells,&Vm[offset], &(s_[offset]),offset,gates_);
   if (id ==3) TT06Func::updateGate3(dt, nCells,&Vm[offset], &(s_[offset]),offset,gates_);
#endif
   profileStop(gateTimer);
}


void TT06Dev_Reaction::initializeMembraneVoltage(std::vector<double>& Vm)
{
   assert(Vm.size() >= s_.size());
   for (unsigned ii=0; ii<s_.size(); ++ii)
      Vm[ii] = TT06Func::defaultVoltage(s_[ii].cellType);
}
