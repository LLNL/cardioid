#ifndef SIMULATE_HH
#define SIMULATE_HH

#include <map>
#include <string>
#include <set>
#include <vector>
#include <cmath>

#include "Long64.hh"
#include "Anatomy.hh"
#include "ThreadServer.hh"
#include "VectorDouble32.hh"

class Diffusion;
class Reaction;
class Stimulus;
class Sensor;
class CommTable;
using std::isnan;

// storage class for persistent data such as potentials that may 
// need to be analyzed or printed out by sensors
struct CheckRange { bool on; double vMin,vMax;};
class PotentialData
{
 public:
 
   PotentialData()
   {
   };
   
   void setup(const Anatomy& anatomy)
   {
      VmArray_     = VectorDouble32(anatomy.size(), 0.);
      dVmDiffusion_= VectorDouble32(anatomy.nLocal(), 0.);
      dVmReaction_ = VectorDouble32(anatomy.nLocal(), 0.);
      unsigned paddedSize = 4*((anatomy.nLocal()+3)/4);
//#include "slow_fix.hh"
//#ifdef SLOW_FIX
      {
	int nFourVecs = paddedSize>>2;
	if(0) paddedSize += 4*((10 - (nFourVecs % 8)) % 8);
	else  paddedSize += ((10 - (nFourVecs & 7)) & 7) << 2;
	VmArray_.reserve(paddedSize);
      }
//#endif
      dVmDiffusion_.reserve(paddedSize);
      dVmReaction_.reserve(paddedSize);

      assert((size_t)&(VmArray_[0])      % 32 == 0);
      assert((size_t)&(dVmDiffusion_[0]) % 32 == 0);
      assert((size_t)&(dVmReaction_[0])  % 32 == 0);
      
   }
   
   // use pointers to vector so that they can be swapped
   VectorDouble32 VmArray_; // local and remote
   VectorDouble32 dVmDiffusion_;
   VectorDouble32 dVmReaction_;
};

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

   enum LoopType {omp, pdr, allSkate, lag};
   
   void checkRanges(int begin, int end,
                    const VectorDouble32& dVmReaction,
                    const VectorDouble32& dVmDiffusion);
   void checkRanges(const VectorDouble32& dVmReaction,
                    const VectorDouble32& dVmDiffusion);
   

   bool checkIO(int loop=-1)const;
   void bufferReactionData(const int begin, const int end);
   void bufferReactionData();
   
   CheckRange checkRange_;
   LoopType loopType_;
   volatile int loop_; // volatile (read + modified in threaded section)
   int maxLoop_;
   int globalSyncRate_;
   int printRate_;
   int64_t printGid_;
   int printRank_;
   int printIndex_;
   FILE *printFile_; 
   int snapshotRate_;
   std::set<Long64> snapshotCellList_;
   int checkpointRate_;
   int nFiles_;

   ThreadTeam diffusionThreads_;
   ThreadTeam reactionThreads_;
   
   double dt_;
   double time_;
   
   int nx_, ny_, nz_; // global size of grid

   std::string name_;
   std::vector<std::string> stateFilename_;

   std::vector<int> sendMap_;
   CommTable* commTable_;
   
   PotentialData vdata_;

   Anatomy anatomy_;
   Diffusion* diffusion_;
   Reaction* reaction_; 
   std::vector<Stimulus*> stimulus_;
   std::vector<Sensor*> sensor_;
   
   void initSensors(const std::vector<std::string>& names);

 private:
   void outOfRange(unsigned index, double dVmr, double dVmd);

};

#endif
