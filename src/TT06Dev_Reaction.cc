#include "TT06Dev_Reaction.hh"
#include <cmath>
#include  <omp.h>
#include "Anatomy.hh"
#include "TT06Func.hh"
#include "pade.hh"
#include "PerformanceTimers.hh"

using namespace std;
using namespace TT06Func;
using namespace PerformanceTimers;

static FILE *tfile[64]; 

TT06Dev_Reaction::TT06Dev_Reaction(const Anatomy& anatomy,double tolerance,int mod, coreGroup *group)
: nCells_(anatomy.nLocal())
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
   group_ = group; 
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
   for (unsigned ii=0; ii<nCells_; ++ii)
   {
      assert(anatomy.cellType(ii) >= 0 && anatomy.cellType(ii) < 256);
      int cellType = ttType_[anatomy.cellType(ii)];
      TT06Func::initState(&(s_[ii]),cellType); 
   }
   for (int i=0;i<omp_get_max_threads();i++) 
   {
        char filename[512]; 
        sprintf(filename,"TT06_%d",i); 
	//tfile[i] =fopen(filename,"w"); 
   }
   
}
TT06Dev_Reaction::~TT06Dev_Reaction()
{
}
void TT06Dev_Reaction::writeStateDev(int loop)
{
int  map[] = { dVK_i , Na_i , Ca_i , Xr1_gate , Xr2_gate , Xs_gate , m_gate , h_gate , j_gate , Ca_ss , d_gate , f_gate , f2_gate , fCass_gate , s_gate , r_gate , Ca_SR , R_prime , jL_gate};
	for (int i=0;i<nStateVar;i++) 
	{
		int k = map[i]; 
		double state = s_[0].state[k];
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
    int coreID,hwThreadID,threadID,nCores,nHwThreads,nThreads;
    groupInfo(group_,coreID,hwThreadID,threadID,nCores,nHwThreads,nThreads); 
    offset=0; 
    int nCellsCore= partition(coreID,nCells_,nCores,offset); 
    int nCells  = partition(hwThreadID,nCellsCore,nHwThreads,offset); 
    return nCells; 
}

void TT06Dev_Reaction::calc(double dt, const vector<double>& Vm, const vector<double>& iStim, vector<double>& dVm)
{
   TT06Func::updateNonGate(dt, nCells_,&Vm[0], &(s_[0]), &dVm[0]);
   TT06Func::updateGate(dt, nCells_,&Vm[0], &(s_[0]));
}
void TT06Dev_Reaction::updateNonGate(double dt, const vector<double>& Vm, vector<double>& dVR)
{
   int offset; 
   profileStart(nonGateTimer);
   int nCells = nonGateWorkPartition(offset); 
   int ompID = omp_get_thread_num(); 
   //fprintf(tfile[ompID],"offset=%d ncells=%d\n",offset,nCells); 
   TT06Func::updateNonGate(dt, nCells,&Vm[offset], &(s_[offset]), &dVR[offset]);
   profileStop(nonGateTimer);
}
void TT06Dev_Reaction::updateGate(double dt, const vector<double>& Vm)
{
   int offset; 
   profileStart(gateTimer);
   int nCells = nonGateWorkPartition(offset); 
   //TT06Func::update_mGate(dt, nCells_,&Vm[0], &(s_[0]));
   //TT06Func::update_hGate(dt, nCells_,&Vm[0], &(s_[0]));
   //TT06Func::update_jGate(dt, nCells_,&Vm[0], &(s_[0]));
   TT06Func::updateGate(dt, nCells,&Vm[offset], &(s_[offset]));
   //TT06Func::update_sGate(dt, nCells_,&Vm[0], &(s_[0]));
   profileStop(gateTimer);
}


void TT06Dev_Reaction::initializeMembraneVoltage(std::vector<double>& Vm)
{
   assert(Vm.size() >= s_.size());
   for (unsigned ii=0; ii<s_.size(); ++ii)
      Vm[ii] = TT06Func::defaultVoltage(s_[ii].cellType);
}
