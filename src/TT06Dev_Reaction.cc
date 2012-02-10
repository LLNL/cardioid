#include "TT06Dev_Reaction.hh"
#include <cmath>
#include "Anatomy.hh"
#include "TT06Func.hh"
#include "pade.hh"

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
      int maxTerms=64; 
      PADE **fit=makeFit(tolerance,V0,V1,deltaV,mod); 
      for (int i=0;fit[i]!=NULL;i++) padeCalc(fit[i],maxTerms,maxTerms,maxCost); 
      writeFit(fit); 
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
      initState(&(s_[ii]),cellType); 
   }
   
}
TT06Dev_Reaction::~TT06Dev_Reaction()
{
}

void TT06Dev_Reaction::calc(double dt, const vector<double>& Vm, const vector<double>& iStim, vector<double>& dVm)
{
   updateNonGate(dt, nCells_,&Vm[0], &(s_[0]), &dVm[0]);
   updateGate(dt, nCells_,&Vm[0], &(s_[0]));
   double c9 = get_c9(); 
   for (unsigned ii=0; ii<nCells_; ++ii) s_[ii].state[K_i] += dt*iStim[ii]*c9 ;
}


void TT06Dev_Reaction::initializeMembraneVoltage(std::vector<double>& Vm)
{
   assert(Vm.size() >= s_.size());
   for (unsigned ii=0; ii<s_.size(); ++ii)
      Vm[ii] = defaultVoltage(s_[ii].cellType);
}
