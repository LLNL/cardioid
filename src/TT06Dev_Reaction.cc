#include "TT06Dev_Reaction.hh"
#include <cmath>
#include "Anatomy.hh"
#include "TT06Func.hh"

using namespace std;


TT06Dev_Reaction::TT06Dev_Reaction(const Anatomy& anatomy,double tolerance,int mod)
: nCells_(anatomy.nLocal())
{
   static bool initialized = false;
   if (! initialized)
   {
      initialized = true;
      initCnst();
   }
   {
      double V0 = -100.0; 
      double V1 =  50.0; 
      double deltaV = 0.1; 
      int maxCost=128; 
      int maxOrder=64; 
      makeFit(tolerance,V0,V1,deltaV,maxOrder,maxCost,mod); 
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
   for (unsigned ii=0; ii<nCells_; ++ii)
   {
      assert(anatomy.cellType(ii) >= 0 && anatomy.cellType(ii) < 256);
      int cellType = ttType_[anatomy.cellType(ii)];
      s_[ii].cellType = cellType; 
      initState(s_[ii].state,cellType); 
   }
   
}
TT06Dev_Reaction::~TT06Dev_Reaction()
{
}

void TT06Dev_Reaction::calc(double dt, const vector<double>& Vm, const vector<double>& iStim, vector<double>& dVm)
{
   assert(nCells_ == dVm.size());


   double c9 = get_c9(); 
   double rates[nStateVar];
   for (unsigned ii=0; ii<nCells_; ++ii)
   {
      computeNonGateRates(1,&(Vm[ii]), &(s_[ii]), &dVm[ii],rates);
      computeGateRates(dt, Vm[ii], s_[ii].state, s_[ii].cellType, rates);
      rates[K_i] += iStim[ii]*c9 ;
      for (int i=0;i<nStateVar;i++) s_[ii].state[i] += dt*rates[i]; 
   }
}


void TT06Dev_Reaction::initializeMembraneVoltage(std::vector<double>& Vm)
{
   assert(Vm.size() >= s_.size());
   for (unsigned ii=0; ii<s_.size(); ++ii)
      Vm[ii] = defaultVoltage(s_[ii].cellType);
}
