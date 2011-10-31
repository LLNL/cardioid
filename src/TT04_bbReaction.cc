#include "TT04_bbReaction.hh"
#include "Anatomy.hh"
#include "IBM_tenTusscher04.hh"

using namespace std;

TT04_bbReaction::TT04_bbReaction(const Anatomy& anatomy)
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

   
   IBM_tenTusscher04_endoLUT::TT04LUT_Init();
   IBM_tenTusscher04_midLUT::TT04LUT_Init();
   IBM_tenTusscher04_epiLUT::TT04LUT_Init();
   
   cellModel_.reserve(nCells_);
   for (unsigned ii=0; ii<nCells_; ++ii)
   {
      assert(anatomy.cellType(ii) >= 0 && anatomy.cellType(ii) < 256);
      int ttType = ttType_[anatomy.cellType(ii)];
      cellModel_.push_back(IBM_tenTusscher04(NULL, ttType));
   }
   
}

void TT04_bbReaction::calc(double dt, const vector<double>& Vm, vector<double>& dVm)
{
   assert(nCells_ == dVm.size());
   
   for (unsigned ii=0; ii<nCells_; ++ii)
   {
      dVm[ii] = cellModel_[ii].Calc(dt, Vm[ii]);
   }
}


