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
#include "slow_fix.hh"
#include "lazy_array.hh"

class Diffusion;
class ReactionManager;
class Stimulus;
class Sensor;
class Drug;
class CommTable;
//using std::isnan;

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
      VmTransport_.resize(anatomy.size());
      dVmDiffusionTransport_.resize(anatomy.nLocal());
      dVmReactionTransport_.resize(anatomy.nLocal());

      
      //unsigned paddedSize = convertActualSizeToBufferSize(anatomy.nLocal());
      //VmArray.reserve(paddedSize);
      //dVmDiffusion.reserve(paddedSize);
      //dVmReaction.reserve(paddedSize);

      //assert((size_t)&(VmArray[0])      % 32 == 0);
      //assert((size_t)&(dVmDiffusion[0]) % 32 == 0);
      //assert((size_t)&(dVmReaction[0])  % 32 == 0);
   }
   
   // use pointers to vector so that they can be swapped
   lazy_array<double> VmTransport_; // local and remote
   lazy_array<double> dVmDiffusionTransport_;
   lazy_array<double> dVmReactionTransport_;
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

   enum LoopType {omp, pdr};
   
   void checkRanges(int begin, int end,
                    ro_mgarray_ptr<double> Vm,
                    ro_mgarray_ptr<double> dVmReaction,
                    ro_mgarray_ptr<double> dVmDiffusion);
   void checkRanges(ro_mgarray_ptr<double> Vm,
                    ro_mgarray_ptr<double> dVmReaction,
                    ro_mgarray_ptr<double> dVmDiffusion);
   

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
   int checkpointRate_;
   bool asciiCheckpoints_;

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
   ReactionManager* reaction_;
   std::vector<Stimulus*> stimulus_;
   std::vector<Sensor*> sensor_;
    
   void initSensors(const std::vector<std::string>& names);

 private:
   void outOfRange(unsigned index, const double Vm, const double dVmr, const double dVmd);

};

#endif
