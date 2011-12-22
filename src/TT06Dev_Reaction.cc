#include "TT06Dev_Reaction.hh"
#include <cmath>
#include "Anatomy.hh"
//#include "TT06_CellML.hh"
#include "TT06DevFit.hh"

using namespace std;

struct TT06DevState
{
   double state[nStateVar];
};




TT06Dev_Reaction::TT06Dev_Reaction(const Anatomy& anatomy)
: nCells_(anatomy.nLocal())
{
   static bool initialized = false;
   if (! initialized)
   {
      initialized = true;
      initCnst();
      Approx(1301,32,0.005); 
   }
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
      for (unsigned jj=0; jj<nStateVar; ++jj) { initState(s_[ii].state,cellType); }
   }
   
}

TT06Dev_Reaction::~TT06Dev_Reaction()
{
}

void TT06Dev_Reaction::calc(double dt, const vector<double>& Vm, const vector<double>& iStim, vector<double>& dVm)
{
   assert(nCells_ == dVm.size());


   int cellType =0; 
   double c9 = get_c9(); 
   for (unsigned ii=0; ii<nCells_; ++ii)
   {
      double dVdt = computeUpdates(dt, Vm[ii], s_[ii].state, cellType);
      dVm[ii] = dVdt; 
      s_[ii].state[K_i] += iStim[ii]*c9*dt ;
   }
}


