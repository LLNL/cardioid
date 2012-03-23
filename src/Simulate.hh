#ifndef SIMULATE_HH
#define SIMULATE_HH

#include <map>
#include <string>
#include <vector>

#include "Anatomy.hh"
#include "Threading.hh"

class Diffusion;
class Reaction;
class Stimulus;
class Sensor;
class CommTable;

// Kitchen sink class for heart simulation.  This is probably a great
// example of how *not* to do design, but we need to get something
// running before we can understand what good design would be.

// This class intentionally exposes all the members as public items.
// This is poor encapsulation, but it seems silly at this point to make
// get/set calls for all of the externals while we are still figuring
// out what we are doing.

class Simulate
{
 public:

   int loop_;
   int maxLoop_;
   int printRate_;
   int snapshotRate_;
   int checkpointRate_;
   int parallelDiffusionReaction_;
   CoreGroup* diffusionGroup_; 
   CoreGroup* reactionGroup_; 
   Threading tinfo_; 
   
   double dt_;
   double time_;
   
   int nx_, ny_, nz_; // global size of grid

   std::string name_;
   std::vector<std::string> stateFilename_;

   std::vector<int> sendMap_;
   CommTable* commTable_;
   
   std::vector<double> VmArray_; // local and remote

   Anatomy anatomy_;
   Diffusion* diffusion_;
   Reaction* reaction_; 
   std::vector<Stimulus*> stimulus_;
   std::vector<Sensor*> sensor_;
};

#endif
