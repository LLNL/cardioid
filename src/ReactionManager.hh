#ifndef REACTIONMANAGER_HH
#define REACTIONMANAGER_HH

#include <map>
#include <vector>
#include <string>
#include "VectorDouble32.hh"

class Reaction;

class ReactionManager
{
 public:
   void calc(double dt,
             const VectorDouble32& Vm,
             const std::vector<double>& iStim,
             VectorDouble32& dVm);
   void updateNonGate(double dt, const VectorDouble32& Vm, VectorDouble32& dVR);
   void updateGate   (double dt, const VectorDouble32& Vm);

   /** Populates the Vm array with some sensible default initial
    * membrane voltage.  Vm will be the parallel to the local cells in
    * the anatomy that was used to create the concrete reaction class. */
   void initializeMembraneState(VectorDouble32& Vm);

   void scaleCurrents(std::vector<double>);  

   void addReaction(const std::string& reactionName);
   void create(const double dt, const Anatomy& anatomy, const ThreadTeam &group, const std::vector<std::string>& scaleCurrents);

 private:
   std::vector<Reaction*> reactions_;
   std::vector<int> reactionIndexFromLocal_;
   std::vector<int> reactionOffsetFromLocal_;
   std::map<int,int> reactionIndexFromAnatomy_;

   std::vector<VectorDouble32> VmPerReaction_;
   std::vector<VectorDouble32> iStimPerReaction_;
   std::vector<VectorDouble32> dVmPerReaction_;
};

#endif
