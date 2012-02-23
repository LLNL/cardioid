#ifndef STATE_LOADER_HH
#define STATE_LOADER_HH

#include <string>

class BucketOfBits;
class Anatomy;

/** Caller must delete returned pointer */
BucketOfBits* loadAndDistributeState(const std::string& filename,
                                     const Anatomy& anatomy);


#endif
