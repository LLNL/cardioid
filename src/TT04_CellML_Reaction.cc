#include "TT04_CellML_Reaction.hh"
#include <cmath>
#include "object_cc.hh"
#include "Anatomy.hh"
#include "reactionFactory.hh"
#include "TT04_CellML.hh"
#include "TT04_CellML_Endo.hh"
#include "TT04_CellML_Mid.hh"
#include "TT04_CellML_Epi.hh"

using namespace std;

struct TT04_CellMLState
{
   double state[17];
};




TT04_CellML_Reaction::TT04_CellML_Reaction(const int numPoints,
                                           const int ttType,
                                           IntegratorType integrator)
: nCells_(numPoints),
          integrator_(integrator)
{  
   cellModel_.reserve(nCells_);
   for (unsigned ii=0; ii<nCells_; ++ii)
   {
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
                                const VectorDouble32& Vm,
                                const vector<double>& iStim,
                                VectorDouble32& dVm)
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
void TT04_CellML_Reaction::initializeMembraneVoltage(VectorDouble32& Vm)
{
   assert(Vm.size() >= s_.size());
   for (unsigned ii=0; ii<s_.size(); ++ii)
      Vm[ii] = cellModel_[ii]->defaultState(0);
}


   
void TT04_CellML_Reaction::forwardEulerIntegrator(
   double dt,
   const VectorDouble32& Vm,
   const vector<double>& iStim,
   VectorDouble32& dVm)
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
   const VectorDouble32& Vm,
   const vector<double>& iStim,
   VectorDouble32& dVm)
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


REACTION_FACTORY(TT04_CellML)(OBJECT* obj, const double, const int numPoints, const ThreadTeam&)
{
   TT04_CellML_Reaction::IntegratorType integrator;
   string tmp;
   objectGet(obj, "integrator", tmp, "rushLarsen");
   if      (tmp == "rushLarsen")   integrator = TT04_CellML_Reaction::rushLarsen;
   else if (tmp == "rushLarson")   integrator = TT04_CellML_Reaction::rushLarsen;
   else if (tmp == "forwardEuler") integrator = TT04_CellML_Reaction::forwardEuler;
   else    assert(false);
   int ttType;
   objectGet(obj, "ttType", ttType, "0");
      
   return new TT04_CellML_Reaction(numPoints, ttType, integrator);
}
