#include "TT06_RRG_Reaction.hh"
#include <cmath>
#include <string>
#include <map>
#include "Anatomy.hh"
#include "TT06_RRG.hh"
#include "BucketOfBits.hh"
#include "units.h"

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

void TT06_RRG_Reaction::initializeMembraneVoltage(std::vector<double>& Vm)
{
   assert(Vm.size() >= cells_.size());
   for (unsigned ii=0; ii<cells_.size(); ++ii)
      Vm[ii] = cells_[ii].defaultVoltage();
}

void TT06_RRG_Reaction::getCheckpointInfo(vector<string>& fieldNames,
                                          vector<string>& fieldUnits) const
{
   TT06_RRG::getCheckpointInfo(fieldNames, fieldUnits);
}

int TT06_RRG_Reaction::getVarHandle(const string& varName) const
{
   return TT06_RRG::getVarHandle(varName);
}

void TT06_RRG_Reaction::setValue(int iCell, int varHandle, double value)
{
   cells_[iCell].setValue(varHandle, value);
}

double TT06_RRG_Reaction::getValue(int iCell, int varHandle) const
{
   return cells_[iCell].getValue(varHandle);
}

void TT06_RRG_Reaction::getValue(int iCell,
                                 const vector<int>& handle,
                                 vector<double>& value) const
{
   cells_[iCell].getValue(handle, value);
}

const string TT06_RRG_Reaction::getUnit(const string& varName) const
{
   return TT06_RRG::getUnit(varName);
}

