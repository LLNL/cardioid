#include "TT06_RRG_Reaction.hh"
#include <cmath>
#include "Anatomy.hh"
#include "TT06_RRG.hh"

using namespace std;

TT06_RRG_Reaction::TT06_RRG_Reaction(const Anatomy& anatomy)
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

   cells_.reserve(anatomy.nLocal());
   for (unsigned ii=0; ii<anatomy.nLocal(); ++ii)
   {
      assert(anatomy.cellType(ii) >= 0 && anatomy.cellType(ii) < 256);
      int ttType = ttType_[anatomy.cellType(ii)];
      cells_.push_back(TT06_RRG(ttType));
   }
}

void TT06_RRG_Reaction::calc(double dt,
                             const vector<double>& Vm,
                             const vector<double>& iStim,
                             vector<double>& dVm)
{
   assert(cells_.size() == dVm.size());

   int nCells = cells_.size();
#pragma omp parallel for
   for (int ii=0; ii<nCells; ++ii)
      dVm[ii] = cells_[ii].calc(dt, Vm[ii], iStim[ii]);
}
