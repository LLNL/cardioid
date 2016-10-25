#include "TT06_RRG_Reaction.hh"
#include <cmath>
#include <string>
#include <map>
#include "Anatomy.hh"
#include "TT06_RRG.hh"
#include "BucketOfBits.hh"
#include "units.h"

using namespace std;

TT06_RRG_Reaction::TT06_RRG_Reaction(const int numPoints, int ttType, const TT06_RRG_ReactionParms& parms)
{
   initConst(parms); 

   cells_.reserve(numPoints);
   for (unsigned ii=0; ii<numPoints; ++ii)
   {
      if (ii == 0) cells_.push_back(TT06_RRG(ttType, &constants_[0]));
      else cells_.push_back(TT06_RRG(ttType));
   }
}

void TT06_RRG_Reaction::calc(double dt,
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
void TT06_RRG_Reaction::initConst(const TT06_RRG_ReactionParms &parms)
{
   constants_[0] = 8314.472;
   constants_[1] = 310;
   constants_[2] = 96485.3415;
   constants_[3] = 0.185;
   constants_[4] = 0.016404;
   //unused      constants_[5] = 10;
   //unused      constants_[6] = 1000;
   //unused      constants_[7] = 1;
   //unused      constants_[8] = 52;
   constants_[9] = 0.03;
   constants_[10] = 5.4;
   constants_[11] = 140;
   constants_[12] = 2;
   constants_[13] = 5.405;
   constants_[14] = 0.153;
   // unused      constants_[15] = 0.098;
   constants_[16] = 14.838;
   constants_[17] = 0.00029;
   constants_[18] = 0.0000398;
   constants_[19] = 0.000592;
   //unused      constants_[20] = 0.294;
   //unused      constants_[21] = 2.724;
   constants_[22] = 1;
   constants_[23] = 40;
   constants_[24] = 1000;
   constants_[25] = 0.1;
   constants_[26] = 2.5;
   constants_[27] = 0.35;
   constants_[28] = 1.38;
   constants_[29] = 87.5;
   constants_[30] = 0.1238;
   constants_[31] = 0.0005;
   constants_[32] = 0.0146;
   constants_[33] = 0.15;
   constants_[34] = 0.045;
   constants_[35] = 0.06;
   constants_[36] = 0.005;
   constants_[37] = 1.5;
   constants_[38] = 2.5;
   constants_[39] = 1;
   constants_[40] = 0.102;
   constants_[41] = 0.0038;
   constants_[42] = 0.00025;
   constants_[43] = 0.00036;
   constants_[44] = 0.006375;
   constants_[45] = 0.2;
   constants_[46] = 0.001;
   constants_[47] = 10;
   constants_[48] = 0.3;
   constants_[49] = 0.4;
   constants_[50] = 0.00025;
   constants_[51] = 0.001094;
   constants_[52] = 0.00005468;
   if (parms.Ko  != 1.0)  constants_[10] = parms.Ko; 
}


void TT06_RRG_Reaction::initializeMembraneVoltage(VectorDouble32& Vm)
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

