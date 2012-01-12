#include "TT04_CellML_Reaction.hh"
#include <cmath>
#include "Anatomy.hh"
#include "TT04_CellML.hh"
#include "TT04_CellML_Endo.hh"
#include "TT04_CellML_Mid.hh"
#include "TT04_CellML_Epi.hh"

using namespace std;

struct TT04_CellMLState
{
   double state[17];
};




TT04_CellML_Reaction::TT04_CellML_Reaction(const Anatomy& anatomy,
                                           IntegratorType integrator)
: nCells_(anatomy.nLocal()),
          integrator_(integrator)
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

   s_.resize(nCells_);
   for (unsigned ii=0; ii<nCells_; ++ii)
   {
      for (unsigned jj=0; jj<17; ++jj)
      {
         s_[ii].state[jj] = cellModel_[ii]->defaultState(jj);
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

   switch (integrator_)
   {
     case forwardEuler:
      forwardEulerIntegrator(dt, Vm, iStim, dVm);
      break;
     case rushLarsen:
      rushLarsenIntegrator(dt, Vm, iStim, dVm);
      break;
     default:
      assert(false);
   }
}
void TT04_CellML_Reaction::initializeMembraneVoltage(std::vector<double>& Vm)
{
   assert(Vm.size() >= s_.size());
   for (unsigned ii=0; ii<s_.size(); ++ii)
      Vm[ii] = cellModel_[ii]->defaultState(0);
}


   
void TT04_CellML_Reaction::forwardEulerIntegrator(
   double dt,
   const vector<double>& Vm,
   const vector<double>& iStim,
   vector<double>& dVm)
{
   for (unsigned ii=0; ii<nCells_; ++ii)
   {
      double rates[17];
      double algebraic[67];
      dVm[ii] = cellModel_[ii]->calc(Vm[ii], iStim[ii], s_[ii].state, rates, algebraic);

      // forward euler to integrate internal state variables.
      for (unsigned jj=1; jj<17; ++jj)
         s_[ii].state[jj] += rates[jj] * dt;
   
   }
}

void TT04_CellML_Reaction::rushLarsenIntegrator(
   double dt,
   const vector<double>& Vm,
   const vector<double>& iStim,
   vector<double>& dVm)
{
   for (unsigned ii=0; ii<nCells_; ++ii)
   {
      double rates[17];
      double algebraic[67];
      dVm[ii] = cellModel_[ii]->calc(Vm[ii], iStim[ii], s_[ii].state, rates, algebraic);

      // forward euler for all states except rushLarsen for fast sodium m gate.
      for (unsigned jj=1; jj<7; ++jj)
         s_[ii].state[jj] += rates[jj] * dt;
      s_[ii].state[7] =  algebraic[4] - (algebraic[4]-s_[ii].state[7])*exp(-dt/algebraic[39]);
      for (unsigned jj=8; jj<17; ++jj)
         s_[ii].state[jj] += rates[jj] * dt;
   
   }
}

