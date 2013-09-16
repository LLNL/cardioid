#include "Reaction.hh"
#include <cassert>

using namespace std;

/** Gets the names and units of all fields that should be written in a
 *  checkpoint file. */
void Reaction::getCheckpointInfo(vector<string>& fieldNames,
                                 vector<string>& fieldUnits) const
{
   // Concrete classes must provide override to enable
   // checkpoint/restart.
   assert(false);
}

/** getVarHandle must return -1 for any name that is not defined in
 *  the concrete reaction class. */
int Reaction::getVarHandle(const string& varName) const
{
   assert(false);
   return 0;
}

/** getVarHandle must return -1 for any name that is not defined in
 *  the concrete reaction class. */
vector<int> Reaction::getVarHandle(const vector<string>& varName) const
{
   vector<int> handle;
   for (unsigned ii=0; ii<varName.size(); ++ii)
      handle.push_back(getVarHandle(varName[ii]));

   return handle;
}

void Reaction::setValue(int iCell, int varHandle, double value)
{
   assert(false);
}


double Reaction::getValue(int iCell, int handle) const
{
   assert(false);
   return 0.0;
}

void Reaction::getValue(int iCell,
                        const vector<int>& handle,
                        vector<double>& value) const
{
   for (unsigned ii=0; ii<handle.size(); ++ii)
      value[ii] = getValue(iCell, handle[ii]);
}

const string Reaction::getUnit(const string& varName) const
{
   assert(false);
   return string();
}

void Reaction::scaleCurrents(std::vector<double>)
{
   assert(false);
   return;
}
