#include "OHaraRudy_Reaction.hh"
#include <cmath>
#include <string>
#include <map>
#include "Anatomy.hh"
#include "OHaraRudy.hh"
#include "OHaraRudy.h"
#include "BucketOfBits.hh"
#include "units.h"

using namespace std;

OHaraRudy_Reaction::OHaraRudy_Reaction(const Anatomy& anatomy)
{
   ttType_.resize(256, -1); 
   ttType_[30] = ENDO_CELL;
   ttType_[31] = ENDO_CELL;
   ttType_[75] = ENDO_CELL;
   ttType_[76] = M_CELL;
   ttType_[77] = EPI_CELL;
   ttType_[100] = ENDO_CELL;
   ttType_[101] = M_CELL;
   ttType_[102] = EPI_CELL;

   cells_.reserve(anatomy.nLocal());
   for (unsigned ii=0; ii<anatomy.nLocal(); ++ii)
   {
      assert(anatomy.cellType(ii) >= 0 && anatomy.cellType(ii) < 256);
      int ttType = ttType_[anatomy.cellType(ii)];
      cells_.push_back(OHaraRudy(ttType));
   }
}

void OHaraRudy_Reaction::calc(double dt,
                             const VectorDouble32& Vm,
                             const vector<double>& iStim,
                             VectorDouble32& dVm)
{
   assert(cells_.size() == dVm.size());

   int nCells = cells_.size();
#pragma omp parallel for
   for (int ii=0; ii<nCells; ++ii)
      dVm[ii] = cells_[ii].calc(dt, Vm[ii], iStim[ii]);
}

void OHaraRudy_Reaction::initializeMembraneVoltage(VectorDouble32& Vm)
{
   assert(Vm.size() >= cells_.size());
   for (unsigned ii=0; ii<cells_.size(); ++ii)
      Vm[ii] = cells_[ii].defaultVoltage();
}

void OHaraRudy_Reaction::getCheckpointInfo(vector<string>& fieldNames,
                                          vector<string>& fieldUnits) const
{
   OHaraRudy::getCheckpointInfo(fieldNames, fieldUnits);
}

int OHaraRudy_Reaction::getVarHandle(const string& varName) const
{
   return OHaraRudy::getVarHandle(varName);
}

void OHaraRudy_Reaction::setValue(int iCell, int varHandle, double value)
{
   cells_[iCell].setValue(varHandle, value);
}

double OHaraRudy_Reaction::getValue(int iCell, int varHandle) const
{
   return cells_[iCell].getValue(varHandle);
}

void OHaraRudy_Reaction::getValue(int iCell,
                                 const vector<int>& handle,
                                 vector<double>& value) const
{
   cells_[iCell].getValue(handle, value);
}

const string OHaraRudy_Reaction::getUnit(const string& varName) const
{
   return OHaraRudy::getUnit(varName);
}

