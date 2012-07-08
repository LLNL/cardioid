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
class Diffusion;
class Reaction;
class Stimulus;
class Sensor;
class CommTable;
using std::isnan;

// storage class for persistent data such as potentials that may 
// need to be analyzed or printed out by sensors
class PotentialData
{
 public:
 
   PotentialData()
   {
      VmArray_     =NULL;
      dVmDiffusion_=NULL;
      dVmReaction_ =NULL;
   };
   
   ~PotentialData()
   {
      if( VmArray_     !=NULL )delete VmArray_;
      if( dVmDiffusion_!=NULL )delete dVmDiffusion_;
      if( dVmReaction_ !=NULL )delete dVmReaction_;
   };
   
   void setup(const Anatomy& anatomy)
   {
      VmArray_     =new std::vector<double>(anatomy.size(), 0.);
      dVmDiffusion_=new std::vector<double>(anatomy.nLocal(), 0.);
      dVmReaction_ =new std::vector<double>(anatomy.nLocal(), 0.);
   }
   
   std::vector<double>* swapdVmReaction(std::vector<double>* const dVmReaction)
   {
      std::vector<double>* tmp=dVmReaction_;
      dVmReaction_=dVmReaction;
      return tmp;
   }
   
   bool outOfRange(const unsigned ii)const
   {
      const double vMax =   60.;
      const double vMin = -110.;
      return ( (*VmArray_)[ii] > vMax || (*VmArray_)[ii] < vMin || isnan((*VmArray_)[ii]) );
   }
   
   // use pointers to vector so that they can be swapped
   std::vector<double>* VmArray_; // local and remote
   std::vector<double>* dVmDiffusion_;
   std::vector<double>* dVmReaction_;
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

   void checkRanges(const std::vector<double>& dVmReaction,
                    const std::vector<double>& dVmDiffusion);

   bool checkIO()const;
   
   bool checkRanges_;
   volatile int loop_; // volatile (read + modified in threaded section)
   int maxLoop_;
   int printRate_;
   long long int printGid_;
   int printRank_;
   int printIndex_;
   int printInit_; 
   FILE *printFile_; 
   int snapshotRate_;
   std::set<Long64> snapshotCellList_;
   int checkpointRate_;
   int parallelDiffusionReaction_;

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
