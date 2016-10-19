
#include <set>
#include <algorithm>
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

class SortByRidxThenAnatomyThenGid {
 public:
   const map<int, int>& ridxFromTag_;
   SortByRidxThenAnatomyThenGid(map<int, int>& ridxFromTag) : ridxFromTag_(ridxFromTag)
   {}
   bool operator()(const AnatomyCell& a, const AnatomyCell& b)
   {
      if (ridxFromTag_.find(a.cellType_) != ridxFromTag_.find(b.cellType_))
      {
         if (ridxFromTag_.find(a.cellType_) == ridxFromTag_.end())
         {
            return false;
         }
         else if (ridxFromTag_.find(b.cellType_) == ridxFromTag_.end())
         {
            return true;
         }
         else
         {
            return ridxFromTag_.find(a.cellType_)->second < ridxFromTag_.find(b.cellType_)->second;
         }
      }
      else if (a.cellType_ != b.cellType_)
      {
         return a.cellType_ < b.cellType_;
      }
      else if (a.gid_ != b.gid_)
      {
         return a.gid_ < b.gid_;
      }

      return false;
   }
};

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

   vector<int> typeFromRidx(numReactions);
   {
      vector<int> reactionReordering(numReactions);
      {
         int cursor=0;
         for (int itype=0; itype<numTypes; ++itype)
         {
            for (int ireaction=0; ireaction<numReactions; ++ireaction) //bottleneck when #reactions large
            {
               string method;
               objectGet(objects[ireaction], "method", method, "");
               if (method == methodTypes_[itype])
               {
                  /*
                    note, we use cursor instead of ireaction here
                    because we're about to reorder the object arrays.
                    This ensures that the typeFromRidx array
                    corresponds to the same things as what the
                    objects[] array and objectNameFromRidx_ arrays
                    point to.

                    If we used ireaction instead, we would have to
                    permute the arrays below, and we don't want to do
                    that.
                  */
                  typeFromRidx[cursor] = itype;
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
   //now all the reactions are sorted by their type
   //here are the invariants.
   for (int ireaction=0; ireaction<numReactions; ++ireaction)
   {
      assert(objectNameFromRidx_[ireaction] == objects[ireaction]->name);
      string method;
      objectGet(objects[ireaction], "method", method, "");
      assert(method == methodTypes_[typeFromRidx[ireaction]]);
   }
      
   //find all the anatomy tags that have been set as reaction models
   map<int, int> ridxFromTag;
   for (int ireaction=0; ireaction<numReactions; ++ireaction)
   {
      vector<int> anatomyTags;
      objectGet(objects[ireaction], "anatomyTags", anatomyTags);
      for (int itag=0; itag<anatomyTags.size(); ++itag)
      {
         if (ridxFromTag.find(anatomyTags[itag]) != ridxFromTag.end())
         {
            assert(0 && "Duplicate anatomy tags within the reaction models");
         }
         ridxFromTag[anatomyTags[itag]] = ireaction;
      }
   }

   vector<AnatomyCell>& cellArray(anatomy.cellArray());
   //sort the anatomy in the correct order.
   {
      SortByRidxThenAnatomyThenGid cmpFunc(ridxFromTag);
      sort(cellArray.begin(),cellArray.end(), cmpFunc);
   }
      
   //count how many of each we have.
   vector<int> countFromRidx(numReactions);
   for (int ireaction=0; ireaction<numReactions; ++ireaction)
   {
      countFromRidx[ireaction] = 0;
   }
   for (int icell=0; icell<cellArray.size(); ++icell)
   {
      const AnatomyCell& cell(cellArray[icell]);
      if (ridxFromTag.find(cell.cellType_) != ridxFromTag.end())
      {
         countFromRidx[ridxFromTag[cell.cellType_]]++;
      }
   }

   //create the reaction objects
   reactions_.resize(numReactions);
   extents_.resize(numReactions+1);
   extents_[0] = 0;
   for (int ireaction=0; ireaction<numReactions; ++ireaction)
   {
      int localSize = countFromRidx[ireaction];
      //reactions_[ireaction] = reactionFactory(objectNameFromRidx_[ii], dt, localSize, group, scaleCurrents);
      extents_[ireaction+1] = extents_[ireaction]+localSize;
      VmPerReaction_.resize(localSize);
      iStimPerReaction_.resize(localSize);
      dVmPerReaction_.resize(localSize);
   }

   //Ok, now we've created the reaction objects.  Now we need to
   //figure out the state variables and the mapping therein.
   
   
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
