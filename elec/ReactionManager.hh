#ifndef REACTIONMANAGER_HH
#define REACTIONMANAGER_HH

#include <map>
#include <vector>
#include <set>
#include <string>
#include "ThreadUtils.hh"
#include "lazy_array.hh"

class Reaction;

class ReactionManager
{
 public:
   void calc(double dt,
             ro_mgarray_ptr<double> Vm,
             ro_mgarray_ptr<double> iStim,
             wo_mgarray_ptr<double> dVm);
   void updateNonGate(double dt, ro_mgarray_ptr<double> Vm, wo_mgarray_ptr<double> dVR);
   void updateGate   (double dt, ro_mgarray_ptr<double> Vm);
   std::string stateDescription() const;

   /** Populates the Vm array with some sensible default initial
    * membrane voltage.  Vm will be the parallel to the local cells in
    * the anatomy that was used to create the concrete reaction class. */
   void initializeMembraneState(wo_mgarray_ptr<double> Vm);

   void addReaction(const std::string& reactionName);
   void create(const double dt, ro_array_ptr<int> cellTypes, const ThreadTeam &group);

   /** Functions needed for checkpoint/restart */
   void getCheckpointInfo(std::vector<std::string>& fieldNames,
                          std::vector<std::string>& fieldUnits) const;
   int getVarHandle(const std::string& varName) const;
   std::vector<int> getVarHandle(const std::vector<std::string>& varName) const;
   void setValue(int iCell, int varHandle, double value);
   double getValue(int iCell, int varHandle) const;
   void getValue(int iCell,
                 ro_array_ptr<int> handle,
                 wo_array_ptr<double> value) const;
   const std::string getUnit(const std::string& varName) const;
   std::vector<int> allCellTypes() const;
   
 private:
   std::vector<std::string> objectNameFromRidx_;
   std::vector<Reaction*> reactions_;
   std::vector<int> extents_;
   lazy_array<int> EindexFromIindex_;
   std::vector<int> IindexFromEindex_;
   
   std::vector<std::string> unitFromHandle_;
   std::map<std::string, int> handleFromVarname_;

   int getRidxFromCell(const int iCell) const;
   bool subUsesHandle(const int ridx, const int handle, int& subHandle, double& myUnitFromTheirUnit) const;
   
   std::vector<std::map<int, std::pair<int, double> > > subHandleInfoFromTypeAndHandle_;

   std::vector<int> typeFromRidx_;
   std::vector<std::string> methodNameFromType_;
   std::set<int> allCellTypes_;
};

#endif
