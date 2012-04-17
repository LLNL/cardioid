#include "TT06Dev_Reaction.hh"
#include <cmath>
#include "Anatomy.hh"
#include "TT06Func.hh"
#include "TT06Gates.h"
#include "pade.hh"
#include "PerformanceTimers.hh"
#include "ThreadServer.hh"
#include <cstdlib>
#include <map>
#include <string>
#include <algorithm>

using namespace std;
using namespace PerformanceTimers;

int workBundle(int index, int nItems, int nGroups , int mod, int& offset);
static FILE *tfile[64]; 

class SortByRank
{
 public:
   SortByRank(const vector<int>& typeMap)
   : map_(typeMap){};
   bool operator()(const AnatomyCell& a, const AnatomyCell& b)
   {
      const int aRank = map_[a.cellType_];
      const int bRank = map_[b.cellType_];
      assert(aRank >= 0 && bRank >= 0);
      if (aRank < bRank)
         return true;
      if (aRank == bRank)
      {
         if (a.cellType_ < b.cellType_)
            return true;
         if (a.cellType_ == b.cellType_)
            return (a.gid_ < b.gid_);
      }
      return false;
   }
      
 private:
   vector<int> map_;
};


TT06Dev_Reaction::TT06Dev_Reaction( Anatomy& anatomy,  map<string,CellTypeParmsFull> cellTypeParmsMap, vector<string> cellTypeNames, double tolerance, int mod, int fastReaction, const ThreadTeam& group)
: nCells_(anatomy.nLocal()),
  group_(group)
{
   static bool initialized = false;
   if (! initialized)
   {
      initialized = true;
      TT06Func::initCnst();
   }
   nCellTypes_ = cellTypeNames.size(); 
   nCellsOfType_.resize(nCellTypes_,0); 
   vector<AnatomyCell>& cells = anatomy.cellArray();

   vector<int> tissueType2Rank(256, -1);
   ttType_.resize(256, -1); 
   initialVm_.resize(nCellTypes_);
   cellTypeParms_.resize(nCellTypes_); 

   double stateInitial[nCellTypes_][nStateVar];
   for (int i=0;i<nCellTypes_;i++)
   {
      int cellType = i; 
      int cellRank = i; 
      string name = cellTypeNames[i];   
      
      vector<int> indices=cellTypeParmsMap[name].anatomyIndices; 
      for (int j=0;j<indices.size();j++)
      {
          int k=indices[j];
          assert(0<=k && k < 256);
          assert(tissueType2Rank[ k] == -1);
          tissueType2Rank[ k] = cellRank;
          ttType_[k]=cellType; 
      }

      cellTypeParms_[cellType].cellType = cellType;  
      cellTypeParms_[cellType].s_switch = cellTypeParmsMap[name].s_switch;  
      cellTypeParms_[cellType].P_NaK    = cellTypeParmsMap[name].P_NaK;  
      cellTypeParms_[cellType].g_Ks     = cellTypeParmsMap[name].g_Ks ;  
      cellTypeParms_[cellType].g_to     = cellTypeParmsMap[name].g_to ;  
      cellTypeParms_[cellType].g_NaL    = cellTypeParmsMap[name].g_NaL;  
      initialVm_[cellType] = cellTypeParmsMap[name].Vm; 
      int stateCnt=0; 
      map<string,STATE>stateMap = cellTypeParmsMap[name].state; 
      for (map<string,STATE>::iterator it=stateMap.begin();it!=stateMap.end();it++)
      {
        STATE stateDefault = it->second; 
        double value = stateDefault.value; 
        int type = stateDefault.type; 
        int index = stateDefault.index; 
        stateInitial[cellType][index]=value; stateCnt++;
      }
      assert(stateCnt==nStateVar); 
   }
   SortByRank sortFunc(tissueType2Rank);
   sort(cells.begin(), cells.end(), sortFunc);

   int nCellBuffer_ =  4*((nCells_+3)/4); 
   double *stateBuffer_=(double *)malloc(sizeof(double)*nStateVar*nCellBuffer_); 
   cellTypeVector_.resize(nCellBuffer_); 

   state_.resize(nStateVar); 
   for (unsigned jj=0; jj<nStateVar; ++jj) state_[jj]= stateBuffer_+jj*nCellBuffer_; 

   vector<int> nCellsOfType(nCellTypes_,0); 
   double c9=TT06Func::get_c9(); 
   for (unsigned ii=0; ii<nCells_; ++ii)
   {
      assert(anatomy.cellType(ii) >= 0 && anatomy.cellType(ii) < 256);
      int cellType = ttType_[anatomy.cellType(ii)];
      cellTypeVector_[ii] = cellType; 
      for (int j=0;j<nStateVar;j++) state_[j][ii]  = stateInitial[cellType][j]; 
      state_[dVK_i][ii] = state_[K_i][ii]/c9+initialVm_[cellType];
      nCellsOfType[cellType]++; 
   }
   
   {
      fastReaction_=fastReaction; 
      if ( (fabs(tolerance/1e-4 - 1.0) < 1e-12) && mod == 1 && (fastReaction_ != 0) ) fastReaction_ = 1; 
      if (fastReaction_ == -1) fastReaction_=0; 
      initExp(); 
      double V0 = -100.0; 
      double V1 =  50.0; 
      double deltaV = 0.1; 
      int maxCost=128; 
      int maxTerms=64; 
      PADE **fit_=TT06Func::makeFit(tolerance,V0,V1,deltaV,mod); 
      for (int i=0;fit_[i]!=NULL;i++) padeCalc(fit_[i],maxTerms,maxTerms,maxCost); 
      TT06Func::writeFit(fit_); 
      PADE **gatefit_=fit_+gateFitOffset; 
      int size=0; 
      int i=0; 
   }
   int nThreads = group_.nThreads();
   int nSquads = group_.nSquads();
   int squadSize =1; 
   if (nSquads >0) squadSize= nThreads/nSquads;
   int nEq = 12/squadSize;
   gateWork_.resize(nThreads); 
   
   //printf ("nThreads=%d nSquads=%d squadSize=%d\n",nThreads,nSquads,squadSize); 
   //for (int id=0;id<omp_get_max_threads();id++) 
   for (int id=0;id<nThreads;id++) 
   {
   	const ThreadRankInfo& rankInfo = group_.rankInfo(id);
        //int teamRank = rankInfo.teamRank_; 
         int teamRank = id; 
        if (teamRank == -1) continue; 
//        int squadID   =rankInfo.coreRank_; 
//      int squadRank =rankInfo.squadRank_; 
        int offset; 
       int squadID    = id/squadSize; 
       int squadRank  = id%squadSize ;
        int nCell = workBundle(squadID, nCells_, nSquads , 4, offset);
        gateWork_[teamRank].offsetCell =  offset; 
        gateWork_[teamRank].nCell =   nCell; 
        gateWork_[teamRank].offsetEq = nEq*squadRank; 
        gateWork_[teamRank].nEq     =  nEq; 
        //printf("%d %d : %d %d %d %d\n",id,squadRank,nCell,offset,gateWork_[id].nEq,gateWork_[id].offsetEq); 
   }
   
   printf("fastReaction=%d\n",fastReaction_); 

   dtForFit_=0.0; 
}
TT06Dev_Reaction::~TT06Dev_Reaction()
{
}
void TT06Dev_Reaction::writeStateDev(int loop)
{
   int  map[] = { dVK_i , Na_i , Ca_i , Xr1_gate , Xr2_gate , Xs_gate , m_gate , h_gate , j_gate , Ca_ss , d_gate , f_gate , f2_gate , fCass , s_gate , r_gate , Ca_SR , R_prime , jL_gate};
   int  isGate[] = {    0  ,     0, 0 , 1 , 1 , 1 , 1 , 1 , 1 , 0 , 1 , 1 , 1 , 0 , 1 , 1 , 0 , 0 , 1};
   for (int i=0;i<nStateVar;i++) 
   {
      int k = map[i]; 
      double state; 
      state = state_[k][0];
      printf("%d %24.14le\n",i,state); 
   }
}
int workBundle(int index, int nItems, int nGroups , int mod, int& offset)
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
int TT06Dev_Reaction::nonGateWorkPartition(int& offset)
{
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
   WORK work ={ 0,nCells_,0,12}; 
   TT06Func::updateNonGate(dt, &(cellTypeParms_[0]), nCells_, &(cellTypeVector_[0]),const_cast<double *>(&Vm[0]),  0, &state_[0], &dVm[0]);
   if (fastReaction_ == 1) TT06Func::updateGateFast(dt, nCells_, &(cellTypeVector_[0]), const_cast<double *>(&Vm[0]), 0, &state_[gateOffset],work);
   else TT06Func::updateGate(dt, nCells_, &(cellTypeVector_[0]), const_cast<double *>(&Vm[0]), 0, &state_[gateOffset],work);
/*
   char filename[64]; 
   sprintf(filename,"Nstate0_%1d",fastReaction_); 
   static FILE *file=NULL; 
   if (file == NULL)  file = fopen(filename,"w");
   //for (int i=0;i<nCells_;i++)
   int i=3333;
   {
      fprintf(file,"%5d",i); 
      for(int j=0;j<nStateVar;j++)
      {
          fprintf(file," %24.14f",state_[j][i]); 
      }
      fprintf(file,"\n"); 
   }
*/
}
void TT06Dev_Reaction::updateNonGate(double dt, const vector<double>& Vm, vector<double>& dVR)
{
   int offset; 
   startTimer(nonGateTimer);
   int nCells = nonGateWorkPartition(offset); 
   TT06Func::updateNonGate(dt, &cellTypeParms_[0], nCells, &cellTypeVector_[offset], const_cast<double *>(&Vm[offset]),  offset, &state_[0], &dVR[offset]);
   stopTimer(nonGateTimer);
}
void TT06Dev_Reaction::updateGate(double dt, const vector<double>& Vm)
{
   startTimer(gateTimer);

   const ThreadRankInfo& rankInfo = group_.rankInfo();
   int teamRank = rankInfo.teamRank_;
   if (fastReaction_ == 1) TT06Func::updateGateFast(dt, nCells_, &(cellTypeVector_[0]), const_cast<double *>(&Vm[0]), 0, &state_[gateOffset],gateWork_[teamRank]);
   else TT06Func::updateGate(dt, nCells_, &(cellTypeVector_[0]), const_cast<double *>(&Vm[0]), 0, &state_[gateOffset],gateWork_[teamRank]);
   stopTimer(gateTimer);
/*
   char filename[64]; 
   sprintf(filename,"Nstate1_%1d",fastReaction_); 
   static FILE *file=NULL; 
   if (file == NULL)  file = fopen(filename,"w");
   //for (int i=0;i<nCells_;i++)
   int i=3333;
   {
      fprintf(file,"%5d",i); 
      for(int j=0;j<nStateVar;j++)
      {
          fprintf(file," %24.14f",state_[j][i]); 
      }
      fprintf(file,"\n"); 
   }
*/
}


void TT06Dev_Reaction::initializeMembraneVoltage(std::vector<double>& Vm)
{
   assert(Vm.size() >= cellTypeVector_.size());
   for (unsigned ii=0; ii<cellTypeVector_.size(); ++ii) Vm[ii] = initialVm_[cellTypeVector_[ii]];
}
