#include "TT04_CellML_Reaction.hh"
#include "Anatomy.hh"
#include "TT04_CellML.hh"
#include "TT04_CellML_Endo.hh"
#include "TT04_CellML_Mid.hh"
#include "TT04_CellML_Epi.hh"

using namespace std;

TT04_CellML_Reaction::TT04_CellML_Reaction(const Anatomy& anatomy)
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
      switch (ttType)
      {
	case 0:
	 cellModel_.push_back(new TT04_CellML_Endo());
	 break;
	case 1:
	 cellModel_.push_back(new TT04_CellML_Mid());
	 break;
	case 2:
	 cellModel_.push_back(new TT04_CellML_Epi());
	 break;
	default:
	 assert(false);
      }
   }
   
}

TT04_CellML_Reaction::~TT04_CellML_Reaction()
{
   for (unsigned ii=0; ii<cellModel_.size(); ++ii)
   {
      delete cellModel_[ii];
   }
}


void TT04_CellML_Reaction::calc(double dt,
				const vector<double>& Vm,
				const vector<double>& iStim,
				vector<double>& dVm)
{
   assert(nCells_ == dVm.size());
   
   for (unsigned ii=0; ii<nCells_; ++ii)
   {
      dVm[ii] = cellModel_[ii]->calc(dt, Vm[ii], iStim[ii]);
   }
}


