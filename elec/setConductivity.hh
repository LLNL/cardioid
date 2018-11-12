#ifndef SET_CONDUCTIVITY_HH
#define SET_CONDUCTIVITY_HH

#include <string>
#include <vector>
#include <set>

class BucketOfBits;
class Tuple;
class AnatomyCell;

void setConductivity(const std::string& name,
                     const BucketOfBits& data,
                     const Tuple& globalGridSize,
                     std::vector<AnatomyCell>& cell);

#endif
