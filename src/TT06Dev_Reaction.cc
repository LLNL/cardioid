#include "TT06Dev_Reaction.hh"
#include <cmath>
#include "Anatomy.hh"
#include "TT06_CellML.hh"
#include "TT06Dev.hh"

using namespace std;

struct TT06DevState
{
   double state[19];
};




TT06Dev_Reaction::TT06Dev_Reaction(const Anatomy& anatomy)
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
      int cellType = ttType_[anatomy.cellType(ii)];
      cellModel_.push_back(TT06Dev(cellType));
   }

   s_.resize(nCells_);
   for (unsigned ii=0; ii<nCells_; ++ii)
   {
      for (unsigned jj=0; jj<19; ++jj)
      {
         s_[ii].state[jj] = cellModel_[ii].defaultState(jj);
      }
   }
   
}

TT06Dev_Reaction::~TT06Dev_Reaction()
{
}

void TT06Dev_Reaction::calc(double dt, const vector<double>& Vm, const vector<double>& iStim, vector<double>& dVm)
{
   assert(nCells_ == dVm.size());


   for (unsigned ii=0; ii<nCells_; ++ii)
   {
      dVm[ii] = cellModel_[ii].calc(dt,Vm[ii], iStim[ii], s_[ii].state);
   }
}

