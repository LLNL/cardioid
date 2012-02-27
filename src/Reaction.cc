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

/** getVarHandle must return 0 for any name that is not defined in
 *  the concrete reaction class. */
int Reaction::getVarHandle(const string& varName) const
{
   assert(false);
}

/** getVarHandle must return 0 for any name that is not defined in
 *  the concrete reaction class. */
vector<int> Reaction::getVarHandle(const vector<string>& varName) const
{
   assert(false);
}

void Reaction::setValue(int iCell, int varHandle, double value)
{
   assert(false);
}


double Reaction::getValue(int iCell, int handle) const
{
   assert(false);
}

void Reaction::getValue(int iCell,
                        const vector<int>& handle,
                        vector<double>& value) const
{
   for (unsigned ii=0; ii<handle.size(); ++ii)
      value[ii] = getValue(iCell, handle[ii]);
}

const string& Reaction::getUnit(const string& varName) const
{
   assert(false);
}

