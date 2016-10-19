
#include <set>
#include "ReactionManager.hh"
#include "Reaction.hh"
#include "object_cc.hh"


using namespace std;


///pass through routines.
void ReactionManager::calc(double dt,
                           const VectorDouble32& Vm,
                           const std::vector<double>& iStim,
                           VectorDouble32& dVm)
{
   for (int ii=0; ii<reactions_.size(); ++ii)
   {
      copy(&Vm[extents_[ii]], &Vm[extents_[ii+1]], VmPerReaction_[ii].begin());
      copy(&iStim[extents_[ii]], &iStim[extents_[ii+1]], iStimPerReaction_[ii].begin());
      reactions_[ii]->calc(dt, VmPerReaction_[ii], iStimPerReaction_[ii], dVmPerReaction_[ii]);
      copy(dVmPerReaction_[ii].begin(), dVmPerReaction_[ii].end(), &dVm[extents_[ii]]);
   }
}
   
void ReactionManager::updateNonGate(double dt, const VectorDouble32& Vm, VectorDouble32& dVR)
{
   for (int ii=0; ii<reactions_.size(); ++ii)
   {
      copy(&Vm[extents_[ii]], &Vm[extents_[ii+1]], VmPerReaction_[ii].begin());
      reactions_[ii]->updateNonGate(dt, VmPerReaction_[ii], dVmPerReaction_[ii]);
      copy(dVmPerReaction_[ii].begin(), dVmPerReaction_[ii].end(), &dVR[extents_[ii]]);
   }
}
void ReactionManager::updateGate(double dt, const VectorDouble32& Vm)
{
   for (int ii=0; ii<reactions_.size(); ++ii)
   {
      copy(&Vm[extents_[ii]], &Vm[extents_[ii+1]], VmPerReaction_[ii].begin());
      reactions_[ii]->updateGate(dt, VmPerReaction_[ii]);
   }
}

void ReactionManager::scaleCurrents(std::vector<double> arg)
{
   for (int ii=0; ii<reactions_.size(); ++ii)
   {
      reactions_[ii]->scaleCurrents(arg);
   }
}

/** Populates the Vm array with some sensible default initial
 * membrane voltage.  Vm will be the parallel to the local cells in
 * the anatomy that was used to create the concrete reaction class. */
void ReactionManager::initializeMembraneState(VectorDouble32& Vm)
{
   for (int ii=0; ii<reactions_.size(); ++ii)
   {
      Reaction* reaction = reactions_[ii];
      string objectName = objectNameFromRidx_[ii];
      ::initializeMembraneState(reaction, objectName, VmPerReaction_[ii]);
      copy(VmPerReaction_[ii].begin(), VmPerReaction_[ii].end(), &Vm[extents_[ii]]);      
   }
}


void ReactionManager::addReaction(const std::string& rxnObjectName)
{
   objectNameFromRidx_.push_back(rxnObjectName);
}

void ReactionManager::create(const double dt, Anatomy& anatomy, const ThreadTeam &group, const std::vector<std::string>& scaleCurrents)
{
   //construct an array of all the objects
   int numReactions=objectNameFromRidx_.size();
   vector<OBJECT*> objects(numReactions);
   for (int ii=0; ii<numReactions; ++ii)
   {
       objects[ii] = objectFind(objectNameFromRidx_[ii], "REACTION");
   }

   //get all the method types
   int numTypes;
   {
      set<string> methodTypeSet;
      for (int ii=0; ii<numReactions; ++ii)
      {
         string method;
         objectGet(objects[ii], "method", method, "");
         assert(method != "");
         methodTypeSet.insert(method);
      }
      numTypes = methodTypeSet.size();
      methodTypes_.resize(numTypes);
      copy(methodTypeSet.begin(), methodTypeSet.end(), methodTypes_.begin());
   }

   {
      vector<int> reactionReordering(numReactions);
      {
         int cursor=0;
         for (int itype=0; itype<numTypes; ++itype)
         {
            for (int ireaction=0; ireaction<numReactions; ++ireaction)
            {
               string method;
               objectGet(objects[ireaction], "method", method, "");
               if (method == methodTypes_[itype])
               {
                  reactionReordering[ireaction] = cursor++;
               }
            }
         }
         assert(cursor == numReactions);
      }

      vector<string> nameCopy(numReactions);
      vector<OBJECT*> objectCopy(numReactions);
      for (int ireaction=0; ireaction<numReactions; ++ireaction)
      {
         objectCopy[reactionReordering[ireaction]] = objects[ireaction];
         nameCopy[reactionReordering[ireaction]] = objectNameFromRidx_[ireaction];
      }
      objects = objectCopy;
      objectNameFromRidx_ = nameCopy;
   }

   //ok, now all the reactions are sorted by their type
   
   /*
   reactions_.resize(objectNameFromRidx_.size());
   extents_.resize(objectNameFromRidx_.size()+1);
   extents_[0] = 0;
      reactions_[ii] = reactionFactory(objectNameFromRidx_[ii], dt, anatomy, group, scaleCurrents);
      extents_[ii+1] = extents_[ii]+localSize;
   }
   */
}

std::string ReactionManager::stateDescription() const {
   assert(methodTypes_.size() >= 1);
   string retval = methodTypes_[0];
   for (int itype=1; itype<methodTypes_.size(); ++itype) {
      retval += "+";
      retval += methodTypes_[itype];
   }
   return retval;
}

   /** Functions needed for checkpoint/restart */
void ReactionManager::getCheckpointInfo(std::vector<std::string>& fieldNames,
                                        std::vector<std::string>& fieldUnits) const
{
   fieldNames.clear();
   fieldNames.reserve(handleFromVarname_.size());
   for(map<string,int>::const_iterator iter=handleFromVarname_.begin();
       iter != handleFromVarname_.end();
       ++iter)
   {
      fieldNames.push_back(iter->first);
   }
   for (int ii=0; ii<fieldNames.size(); ++ii)
   {
      fieldUnits.push_back(getUnit(fieldNames[ii]));
   }
}
int ReactionManager::getVarHandle(const std::string& varName) const
{
   assert(handleFromVarname_.find(varName) != handleFromVarname_.end());
   return handleFromVarname_.find(varName)->second;
}

vector<int> ReactionManager::getVarHandle(const vector<string>& varName) const
{
   vector<int> handle;
   for (unsigned ii=0; ii<varName.size(); ++ii)
      handle.push_back(getVarHandle(varName[ii]));

   return handle;
}
void ReactionManager::setValue(int iCell, int varHandle, double value)
{
   int ridx = getRidxFromCell(iCell);
   int subHandle;
   double myUnitFromTheirUnit;
   if (subUsesHandle(ridx, varHandle, subHandle, myUnitFromTheirUnit)) {
      int subCell = iCell-extents_[ridx];
      reactions_[ridx]->setValue(subCell, subHandle, value/myUnitFromTheirUnit);
   }
}
double ReactionManager::getValue(int iCell, int varHandle) const
{
   int ridx = getRidxFromCell(iCell);
   int subHandle;
   double myUnitFromTheirUnit;
   if (subUsesHandle(ridx, varHandle, subHandle, myUnitFromTheirUnit)) {
      int subCell = iCell-extents_[ridx];
      return myUnitFromTheirUnit*reactions_[ridx]->getValue(subCell, subHandle);
   } else {
      return numeric_limits<double>::quiet_NaN();
   }
}
void ReactionManager::getValue(int iCell,
                        const vector<int>& handle,
                        vector<double>& value) const
{
   for (unsigned ii=0; ii<handle.size(); ++ii)
      value[ii] = getValue(iCell, handle[ii]);
}
const std::string ReactionManager::getUnit(const std::string& varName) const
{
   return unitFromHandle_[getVarHandle(varName)];
}
   
int ReactionManager::getRidxFromCell(const int iCell) const {
   //binary search to find the proper Ridx.
   int begin=0;
   int end=extents_.size();
   while (begin+1 < end)
   {
      int mid=(begin+end)/2;
      if (extents_[mid] == iCell)
      {
         return mid;
      }
      else if (extents_[mid] < iCell)
      {
         begin=mid;
      }
      else
      {
         end=mid;
      }
   }
   return begin;
}

bool ReactionManager::subUsesHandle(const int ridx, const int handle, int& subHandle, double& myUnitFromTheirUnit) const
{
   if (subHandleInfoFromRidxAndHandle_[ridx].find(handle) == subHandleInfoFromRidxAndHandle_[ridx].end())
   {
      return false;
   }

   subHandle = subHandleInfoFromRidxAndHandle_[ridx].find(handle)->second.first;
   myUnitFromTheirUnit = subHandleInfoFromRidxAndHandle_[ridx].find(handle)->second.second;
   return true;
}
