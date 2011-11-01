#include "TT04Reaction.hh"
#include "Anatomy.hh"
#include "IBM_TT04.hh"

using namespace std;

TT04Reaction::TT04Reaction(const Anatomy& anatomy)
: nCells_(anatomy.nLocal())
{
   ttType_.resize(256, -1); 
   ttType_[30] = 0;
   ttType_[31] = 0;
   ttType_[75] = 0;
   ttType_[76] = 1;
   ttType_[77] = 2;
   ttType_[100] = 0;
   ttType_[101] = 1;
   ttType_[102] = 2;

   
   cellModel_.reserve(nCells_);
   for (unsigned ii=0; ii<nCells_; ++ii)
   {
      assert(anatomy.cellType(ii) >= 0 && anatomy.cellType(ii) < 256);
      int ttType = ttType_[anatomy.cellType(ii)];
      cellModel_.push_back(IBM_TT04(NULL, ttType));
   }
   
}

void TT04Reaction::calc(double dt,
			const vector<double>& Vm,
			const vector<double>& iStim,
			vector<double>& dVm)
{
   assert(nCells_ == dVm.size());
   
   for (unsigned ii=0; ii<nCells_; ++ii)
   {
      dVm[ii] = - cellModel_[ii].Calc(1000*dt, 1000*Vm[ii], iStim[ii]);
   }
}


