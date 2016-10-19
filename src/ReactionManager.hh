#ifndef REACTIONMANAGER_HH
#define REACTIONMANAGER_HH

#include <map>
#include <vector>
#include <string>
#include "VectorDouble32.hh"
#include "ThreadUtils.hh"
#include "Anatomy.hh"

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
   std::string stateDescription() const;

   /** Populates the Vm array with some sensible default initial
    * membrane voltage.  Vm will be the parallel to the local cells in
    * the anatomy that was used to create the concrete reaction class. */
   void initializeMembraneState(VectorDouble32& Vm);

   void scaleCurrents(std::vector<double>);  

   void addReaction(const std::string& reactionName);
   void create(const double dt, Anatomy& anatomy, const ThreadTeam &group, const std::vector<std::string>& scaleCurrents);

   /** Functions needed for checkpoint/restart */
   void getCheckpointInfo(std::vector<std::string>& fieldNames,
                          std::vector<std::string>& fieldUnits) const;
   int getVarHandle(const std::string& varName) const;
   std::vector<int> getVarHandle(const std::vector<std::string>& varName) const;
   void setValue(int iCell, int varHandle, double value);
   double getValue(int iCell, int varHandle) const;
   void getValue(int iCell,
                 const std::vector<int>& handle,
                 std::vector<double>& value) const;
   const std::string getUnit(const std::string& varName) const;
   
 private:
   std::vector<std::string> objectNameFromRidx_;
   std::vector<Reaction*> reactions_;
   std::vector<int> extents_;
   
   std::vector<VectorDouble32> VmPerReaction_;
   std::vector<std::vector<double> > iStimPerReaction_;
   std::vector<VectorDouble32> dVmPerReaction_;

   std::vector<std::string> unitFromHandle_;
   std::map<std::string, int> handleFromVarname_;

   int getRidxFromCell(const int iCell) const;
   bool subUsesHandle(const int ridx, const int handle, int& subHandle, double& myUnitFromTheirUnit) const;
   
   std::vector<std::map<int, std::pair<int, double> > > subHandleInfoFromRidxAndHandle_;
};

#endif
