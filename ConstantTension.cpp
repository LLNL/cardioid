#include "ConstantTension.hpp"
#include <cassert>
#include <cmath>

using namespace std;

namespace ConstantTension
{

string ThisModel::name() const
{
   return "ConstantTension";
}

enum varHandles
{
   _handle_stretch,
   _handle_stretchVel,
   _handle_tension,
   _handle_dtension_dstretchVel,
   _handle_usedTension,
   NUM_HANDLES
};

const string varNames[] = {
   "stretch",
   "stretchVel",
   "tension",
   "dtension_dstretchVel",
   "usedTension",
   ""
};

const string varUnits[] = {
   "1",
   "1/ms",
   "kPa",
   "kPa/ms",
   "kPa",
   ""
};

enum inputOrder
{
   _inIdx_stretch,
   _inIdx_stretchVel,
   NUM_INPUTS
};

enum outputOrder
{
   _outIdx_tension,
   _outIdx_dtension_dstretchVel,
   NUM_OUTPUTS
};

ThisModel::ThisModel(const int numPoints)
{
   _numPoints = numPoints;
}
   
int ThisModel::getHandle(const string& varName) const
{
   for (int ii=0; ii<NUM_HANDLES; ii++)
   {
      if (varName == varNames[ii])
      {
         return ii;
      }
   }
   return -1;
}

string ThisModel::getUnit(const int varHandle) const
{
   assert(varHandle > 0 && varHandle < NUM_HANDLES);
   return varUnits[varHandle];
}

vector<string> ThisModel::getVarGroup(const std::string type) const
{
   vector<string> names;
   if (0) {}
   else if (type == "input")
   {
      names.resize(NUM_INPUTS);
      names[_inIdx_stretch] = "stretch";
      names[_inIdx_stretchVel] = "stretchVel";
   }
   else if (type == "output")
   {
      names.resize(NUM_OUTPUTS);
      names[_outIdx_tension] = "tension";
      names[_outIdx_dtension_dstretchVel] = "dtension_dstretchVel";
   }
   return names;
}


double ThisModel::get(const int varHandle) const
{
   if (0) {}
   else if (varHandle == _handle_usedTension)
   {
      return usedTension;
   }
   return NAN;
}
void ThisModel::set(const int varHandle, const double value)
{
   if (0) {}
   else if (varHandle == _handle_usedTension)
   {
      usedTension = value;
   }
   assert(0 && "Can't set a value for parameter that doesn't exist");
}

double ThisModel::get(const int varHandle, const int iCell) const
{
   if (0) {}
   return NAN;
}

void ThisModel::set(const int varHandle, const int iCell, const double value)
{
   if (0) {}
   assert(0 && "Can't set a value for state variable that doesn't exist");
}

double ThisModel::get(const int varHandle, const int iCell, const double* inputs[]) const
{
   if (0) {}
   return NAN;
}

void ThisModel::resetDefaultParameters()
{
   usedTension = 0;
}

void ThisModel::initialize(const double* const inputs[])
{
}

void ThisModel::tryTimestep(const double dt, const double* const inputs[])
{
}

void ThisModel::outputForNextTimestep(const double dt, const double* const inputs[], double* const outputs[])
{
   for (int icell=0; icell<_numPoints; icell++)
   {
      outputs[_outIdx_tension][icell] = usedTension;
      outputs[_outIdx_dtension_dstretchVel][icell] = 0;
   }
}

void ThisModel::commitTimestep()
{
}
   
}
