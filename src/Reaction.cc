#include "Reaction.hh"
#include <cassert>

using namespace std;

void Reaction::loadState(const BucketOfBits& data)
{
   // loading state data will fail unless the concrete Reaction class
   // provides an override.
   assert(false);
}
