#include "TT06Dev_Reaction.hh"
#include <cmath>
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

#include "slow_fix.hh"

typedef void  (*UPDATEGATE)(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a) ; 

using namespace std;
using namespace PerformanceTimers;
using namespace TT06Func;
int workBundle(int index, int nItems, int nGroups , int mod, int& offset);

static   UPDATEGATE gateEq[] ={ update_mGate_v1, update_hGate_v1, update_jGate_v1, update_Xr1Gate_v1, 
                               update_Xr2Gate_v1, update_XsGate_v1, update_rGate_v1, update_dGate_v1, 
                               update_fGate_v1, update_f2Gate_v1,  update_jLGate_v1, update_s0Gate_v1, 
                               update_s1Gate_v1, update_sGate_v1} ;
static double *mhuX[14]; 
static double *tauRX[14]; 
static double *gateX[14]; 
static UPDATEGATE gateEqX[14]; 
static int gateThreadMapS0[] = { 0,1,2,3,4,5,6,7,8,9,10,11}; 
static int gateThreadMapS1[] = { 0,1,2,3,4,5,6,7,8,9,10,12}; 
static int gateThreadMapS[] =  { 0,1,2,3,4,5,6,7,8,9,10,13}; 


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
      if (aRank <  bRank) return true;
      if (aRank == bRank) return (a.gid_ < b.gid_);
      return false;
   }
      
 private:
   vector<int> map_;
};

class SortByRank_andTuple_factory
{
  /*
    This class has the same major sorting as SortByRank obove,
    but introduces a minor sorting corresponding to the memory
    layout of the voltage block data. This improves the the over
    all simulation loop time by about a microsecond.
  */
  
  struct mytup {
    long long int xc,yc,zc;
    mytup() : xc(0),yc(0),zc(0) {}
    mytup(const Tuple& t) : xc(t.x()),yc(t.y()),zc(t.z()) {}
    long long int x() { return xc; }
    long long int y() { return yc; }
    long long int z() { return zc; }
  };

 public:
  SortByRank_andTuple_factory(const vector<int>& typeMap,const Anatomy& a)
    : map_(typeMap) 
  {

    for(unsigned int i = 0; i<a.size(); i++) 
    {
      tmap[(long long int) a.gid(i)] = mytup(a.globalTuple(i));
    }
  }

  class compare_junk 
  {
    private:
    vector<int> &map_;
    map<long long int,mytup>& tmap;

    public:
    compare_junk(vector<int>& map_in,map<long long int,mytup>& tmapin) : map_(map_in),tmap(tmapin) {}

    bool operator()(const AnatomyCell& a, const AnatomyCell& b)
    {
      const int aRank = map_[a.cellType_];
      const int bRank = map_[b.cellType_];
      assert(aRank >= 0 && bRank >= 0);
      if (aRank <  bRank) return true;
      if (aRank == bRank)
      {
        mytup ta = tmap[a.gid_], tb = tmap[b.gid_];
        if(ta.x() != tb.x()) return ta.x() < tb.x();
        if(ta.y() != tb.y()) return ta.y() < tb.y();
        return ta.z() < tb.z();
      }
      return false;
    }
  };
  
  compare_junk emit_sorter() { return compare_junk(map_,tmap); }

 private:
   vector<int> map_;
   map<long long int,mytup> tmap;
};


// Can't pass parms by const reference since the parms contain a map and
// we access the map with the subscript operator.  This is a non-const
// operation.  Maybe we can re-write this code to avoid the subscript
// operator. 
TT06Dev_Reaction::TT06Dev_Reaction(double dt, Anatomy& anatomy, TT06Dev_ReactionParms& parms, const ThreadTeam& group)
                 : nCells_(anatomy.nLocal()), group_(group)
{
   static bool initialized = false;
   if (! initialized)
   {
      initialized = true;
      TT06Func::initCnst();
      initExp(); 
   }
   nCellBuffer_ =  4*((nCells_+3)/4); 
   int nFourVecs = nCellBuffer_ >> 2;     // Number of full four vectors. 
   if(0) nCellBuffer_ += 4*((10 - (nFourVecs % 8)) % 8);
   else  nCellBuffer_ += ((10 - (nFourVecs & 7)) & 7) << 2;   

   mkCellTypeParms_(anatomy,parms);
   mkCellTypeVector_(anatomy,parms);
   mkState_(parms);
   currentMap_.resize(parms.currentNames.size()); 
   int nCurrents=0; 
   while (currentNames[nCurrents] != "") {nCurrents++;}
   //printf("current sizes: %d %d %d\n",nCurrents,sizeof(double),sizeof(currentScales_)); 
   assert(nCurrents*sizeof(double) == sizeof(currentScales_)); 
   currentScales_=currentScalesDefault; 
   
   for (int i=0;i<parms.currentNames.size();i++)
   {
      int j = 0;
      for (;j<nCurrents;j++) 
         if ( std::string(currentNames[j]) == parms.currentNames[i]) break; 
      assert(j<nCurrents); 
      currentMap_[i]   = j;
   }


   fastReaction_ = 0; 
   if (parms.fastReaction >= 0) ;
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
   for (eq=0;eq<13;eq++) 
   {
      gateEqX[eq]  = gateEq[eq];  
      if (eq < 12) gateX[eq]    = gate[eq];  
      else gateX[12]    = gate[11];  
      mhuX[eq]     = gateFit[2*eq+0].coef; 
      tauRX[eq]    = gateFit[2*eq+1].coef; 
   }
   mkWorkBundles_(parms);
}
void TT06Dev_Reaction::scaleCurrents(vector <double> currentScales)
{
   assert(currentScales.size() == currentMap_.size()); 
   currentScales_=currentScalesDefault; 
   for (int i=0;i<currentScales.size();i++) 
   {
      int j = currentMap_[i];
      ((double *)&currentScales_)[j]=currentScales[i]; 
   }
}
void TT06Dev_Reaction::mkCellTypeParms_(Anatomy& anatomy,TT06Dev_ReactionParms& parms)
{

   nCellTypes_ = parms.cellTypeNames.size(); 
   tissueType2Rank_.resize(256,-1); 

   ttType_.resize(256, -1); 
   initialVm_.resize(nCellTypes_);
   cellTypeParms_.resize(nCellTypes_); 

   stateInitial_.resize(nCellTypes_);
   for (int i=0;i<nCellTypes_;i++) stateInitial_[i].resize(nStateVar);
   for (int i=0;i<nCellTypes_;i++)
   {
      string name = parms.cellTypeNames[i];   
      int cellType = i; 
      int cellRank = i + nCellTypes_ * parms.cellTypeParms[name].s_switch; 

      vector<int> indices=parms.cellTypeParms[name].anatomyIndices; 
      for (int j=0;j<indices.size();j++)
      {
         int k=indices[j];
         assert(0<=k && k < 256);
         assert(tissueType2Rank_[ k] == -1);
         tissueType2Rank_[ k] = cellRank;
         ttType_[k]=cellType; 
      }

      cellTypeParms_[cellType].cellType = cellType;  
      cellTypeParms_[cellType].s_switch = parms.cellTypeParms[name].s_switch;  
      cellTypeParms_[cellType].P_NaK    = parms.cellTypeParms[name].P_NaK;  
      cellTypeParms_[cellType].g_Ks     = parms.cellTypeParms[name].g_Ks ;  
      cellTypeParms_[cellType].g_Kr     = parms.cellTypeParms[name].g_Kr ;  
      cellTypeParms_[cellType].g_to     = parms.cellTypeParms[name].g_to ;  
      cellTypeParms_[cellType].g_NaL    = parms.cellTypeParms[name].g_NaL;  
      initialVm_[cellType] = parms.cellTypeParms[name].Vm; 
      int stateCnt=0; 
      map<string,STATE>stateMap = parms.cellTypeParms[name].state; 
      for (map<string,STATE>::iterator it=stateMap.begin();it!=stateMap.end();it++)
      {
         STATE stateDefault = it->second; 
         double value = stateDefault.value; 
         int type = stateDefault.type; 
         int index = stateDefault.index; 
         stateInitial_[cellType][index]=value; stateCnt++;
      }
      assert(stateCnt==nStateVar); 
   }
}
void TT06Dev_Reaction::mkCellTypeVector_(Anatomy& anatomy, TT06Dev_ReactionParms& parms)
{
   cellTypeVector_.reserve(nCellBuffer_); 
   cellTypeVector_.resize(nCells_); 
   nCellsOfType_.resize(nCellTypes_,0); 
   vector<AnatomyCell>& cells = anatomy.cellArray();

   SortByRank_andTuple_factory sortFunc(tissueType2Rank_,anatomy);
   sort(cells.begin(), cells.end(), sortFunc.emit_sorter());

   nCell_s0_=0; 
   for (unsigned ii=0; ii<nCells_; ++ii) 
   {
      assert(anatomy.cellType(ii) >= 0 && anatomy.cellType(ii) < 256);
      int cellType = ttType_[anatomy.cellType(ii)];
      cellTypeVector_[ii] = cellType; 
      nCellsOfType_[cellType]++; 
      int s0 = cellTypeParms_[cellType].s_switch;
      if (s0 != 0) nCell_s0_++; 
   }

   for (unsigned ii=nCells_; ii<nCellBuffer_;ii++) 
   {
      cellTypeVector_[ii] = cellTypeVector_[ii-1]; 
   }
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
   for (unsigned ii=0; ii<nCells_; ++ii)
   {
      int cellType = cellTypeVector_[ii];
      for (int j=0;j<nStateVar;j++) state_[j][ii]  = stateInitial_[cellType][j]; 
      state_[dVK_i][ii] = state_[K_i][ii]/c9+initialVm_[cellType];
   } 
   for (unsigned ii=nCells_; ii<nCellBuffer_;ii++) 
      for (int j=0;j<nStateVar;j++) state_[j][ii]  =  state_[j][ii-1];  
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
         gateThreadMapS0[ii]  = m; 
         gateThreadMapS1[ii]  = m; 
         gateThreadMapS [ii]  = m; 
         if ( m == 11) 
         {
            gateThreadMapS1[ii]  = 12; 
            gateThreadMapS[ii]   = 13; 
            sGateIndex = ii; 
         }
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
      int *gateThreadMap; 
      int offsetEq = nEq*squadRank; 
      if (offset+nCell >= nCell_s0_ )  gateThreadMap = gateThreadMapS0; 
      else if (  offset <= nCell_s0_ ) gateThreadMap = gateThreadMapS1; 
      else  
      {
         if (  offsetEq <=  sGateIndex && sGateIndex < offsetEq +nEq) 
         {
            assert(splitRank == -1);   //At most  one thread/rank can gave both s0 and s1 cells 
            int nCell0 = nCell_s0_ - offset; 
            int nCell1 = nCell - nCell0; 
            int eq = 12; 
            PADE *gateFit = fit_+gateFitOffset; 
            double *s1Mhu     = gateFit[2*eq+0].coef; 
            double *s1TauR    = gateFit[2*eq+1].coef; 
            sGateInit(s1Mhu, s1TauR, nCell0, nCell1);
            gateThreadMap = gateThreadMapS; 
            splitRank == teamRank; 
         }
      }

      gateWork_[teamRank].offsetCell =  offset; 
      gateWork_[teamRank].nCell    =  nCell; 
      gateWork_[teamRank].nEq      = nEq; 
      gateWork_[teamRank].map      = gateThreadMap+offsetEq; 

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

void TT06Dev_Reaction::calc(double dt, const VectorDouble32& Vm, const vector<double>& iStim, VectorDouble32& dVm)
{
   WORK work ={ 0,nCells_,0,12}; 
   if (nCells_ > 0) 
   {
      update_nonGate_((void*)fit_, &currentScales_,dt,&(cellTypeParms_[0]), nCells_, &(cellTypeVector_[0]),const_cast<double *>(&Vm[0]),  0, &state_[0], &dVm[0]);
      update_gate_(dt, nCells_, &(cellTypeVector_[0]), const_cast<double *>(&Vm[0]), 0, &state_[gateOffset],fit_,work);
   }
}
void TT06Dev_Reaction::updateNonGate(double dt, const VectorDouble32& Vm, VectorDouble32& dVR)
{
   int offset,nCells; 

#ifdef LEGACY_NG_WORKPARTITION
   nCells = nonGateWorkPartition(offset);
#else
   int ompId =  omp_get_thread_num(); 
   const int tid = group_.teamRank();
   offset = nonGateWork_[tid].offsetCell;
   nCells = nonGateWork_[tid].nCell;
#endif

   if(PRINT_WP) {
      static int inited[64] = {0};
      const int tid = group_.teamRank();

      if(inited[tid] == 0) {
         int pid,np;
         const int nt = group_.nThreads();
         static volatile int tag[64] = {0};

         MPI_Comm_size(MPI_COMM_WORLD,&np);
         MPI_Comm_rank(MPI_COMM_WORLD,&pid);	  

         if(tid == 0) {
            int i;

            MPI_Barrier(MPI_COMM_WORLD);
            for(i = 0; i<np; i++) {
               if(i == pid) {
                  int j;
                  printf("@pid=%03d,tid=%02d :NGWP: offset = %6d  nCells = %6d  end = %6d\n",
                        pid,tid,offset,nCells,offset+nCells);

                  tag[tid] = 1;
                  while(tag[nt-1] == 0) {}
                  for(j = 0; j<nt; j++) tag[j] = 0;
               }
               MPI_Barrier(MPI_COMM_WORLD);
            }
         } else {
            while(tag[tid-1] == 0) {}

            printf("@pid=%03d,tid=%02d :NGWP: offset = %6d  nCells = %6d  end = %6d\n",
                  pid,tid,offset,nCells,offset+nCells);

            tag[tid] = 1;
            while(tag[0] > 0) {}
         }

         inited[tid] = 1;
      }
   }
   if (nCells > 0) 
   {
      startTimer(nonGateTimer);
      update_nonGate_(fit_,&currentScales_,dt,&cellTypeParms_[0], nCells, &cellTypeVector_[offset], const_cast<double *>(&Vm[offset]),  offset, &state_[0], &dVR[offset]);
      stopTimer(nonGateTimer);
   }
}
void TT06Dev_Reaction::updateGate(double dt, const VectorDouble32& Vm)
{
   int teamRank = group_.teamRank();

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
      gateEqX[eq](dt, nCell , vm , gateX[eq]+offsetCell, mhuX[eq], tauRX[eq]);
   }
   stopTimer(gateTimer);
}

void TT06Dev_Reaction::initializeMembraneVoltage(VectorDouble32& Vm)
{
   assert(Vm.size() >= cellTypeVector_.size());
   for (unsigned ii=0; ii<cellTypeVector_.size(); ++ii) Vm[ii] = initialVm_[cellTypeVector_[ii]];
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

