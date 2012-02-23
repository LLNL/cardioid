#include "Reaction.hh"
#include <cassert>

using namespace std;

void Reaction::loadState(const BucketOfBits& data)
{
   // loading state data will fail unless the concrete Reaction class
   // provides an override.
   assert(false);
}

/** Gets the names and units of all fields that should be written in a
 *  checkpoint file. */
void Reaction::getCheckpointInfo(vector<string>& fieldNames,
                                 vector<string>& fieldUnits) const
{
   // Concrete classes must provide override to enable
   // checkpoint/restart.
   assert(false);
}

vector<int> Reaction::getHandle(const vector<string>& varName) const
{
   assert(false);
}

void Reaction::getValue(int iCell,
                        const vector<int>& handle,
                        vector<double>& value) const
{
   assert(false);
}
