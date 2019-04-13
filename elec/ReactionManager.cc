
#include <set>
#include <algorithm>
#include "ReactionManager.hh"
#include "Reaction.hh"
#include "object_cc.hh"
#include "units.h"
#include "reactionFactory.hh"
#include "slow_fix.hh"


using namespace std;

///pass through routines.
void ReactionManager::calc(double dt,
                           ro_mgarray_ptr<double> Vm,
                           ro_mgarray_ptr<double> iStim,
                           wo_mgarray_ptr<double> dVm)
{
   for (int ii=0; ii<reactions_.size(); ++ii)
   {
      reactions_[ii]->calc(dt,
                           ro_mgarray_ptr<int>(EindexFromIindex_).slice(extents_[ii],extents_[ii+1]),
                           Vm,
                           iStim,
                           dVm);
   }
}

void ReactionManager::updateNonGate(double dt,
                                    ro_mgarray_ptr<double> Vm,
                                    wo_mgarray_ptr<double> dVm)
{
   for (int ii=0; ii<reactions_.size(); ++ii)
   {
      reactions_[ii]->updateNonGate(dt,
                                    ro_mgarray_ptr<int>(EindexFromIindex_).slice(extents_[ii],extents_[ii+1]),
                                    Vm,
                                    dVm);
   }
}

void ReactionManager::updateGate(double dt,
                                 ro_mgarray_ptr<double> Vm)
{
   for (int ii=0; ii<reactions_.size(); ++ii)
   {
      reactions_[ii]->updateGate(dt,
                                 ro_mgarray_ptr<int>(EindexFromIindex_).slice(extents_[ii],extents_[ii+1]),
                                 Vm);
   }
}

/** Populates the Vm array with some sensible default initial
 * membrane voltage.  Vm will be the parallel to the local cells in
 * the anatomy that was used to create the concrete reaction class. */
void ReactionManager::initializeMembraneState(wo_mgarray_ptr<double> Vm)
{
   for (int ii=0; ii<reactions_.size(); ++ii)
   {
      ::initializeMembraneState(reactions_[ii], objectNameFromRidx_[ii],
                                ro_mgarray_ptr<int>(EindexFromIindex_).slice(extents_[ii],extents_[ii+1]),
                                Vm);
   }
}


void ReactionManager::addReaction(const std::string& rxnObjectName)
{
   objectNameFromRidx_.push_back(rxnObjectName);
}


void ReactionManager::create(const double dt, ro_array_ptr<int> cellTypes, const ThreadTeam &group)
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
      methodNameFromType_.resize(numTypes);
      copy(methodTypeSet.begin(), methodTypeSet.end(), methodNameFromType_.begin());
   }

   typeFromRidx_.resize(numReactions);
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
               if (method == methodNameFromType_[itype])
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
                  typeFromRidx_[cursor] = itype;
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
      assert(method == methodNameFromType_[typeFromRidx_[ireaction]]);
   }

   //check that all the reaction objects we're using have cellTypes defined and that cellTypes are numbers.
   {
      bool foundCellTypeProblem=false;
      for (int ireaction=0; ireaction<numReactions; ++ireaction)
      {
         vector<string> cellTypesText;
         objectGet(objects[ireaction], "cellTypes", cellTypesText);
         if (cellTypesText.empty()) {
            printf("cellTypes not found for %s.  cellTypes now required on all REACTION objects\n",
                   objectNameFromRidx_[ireaction].c_str());
            foundCellTypeProblem=true;
         }
         for (int itag=0; itag<cellTypesText.size(); ++itag)
         {
            for (int cursor=0; cursor<cellTypesText[itag].size(); ++cursor) {
               unsigned char thisChar = cellTypesText[itag][cursor];
               if (!('0' <= thisChar && thisChar <= '9')) {
                  printf("Expected a number, found %s instead for cellTypes attribute on object %s\n"
                         "If you're converting an old object.data file, use cellTypeName instead.\n"
                         "Notice the singular.  You'll need diffent reaction objects for each type.\n",
                         cellTypesText[itag].c_str(), objectNameFromRidx_[ireaction].c_str());
                  foundCellTypeProblem=true;
               }
            }
         }
      }
      assert(!foundCellTypeProblem);
   }
   
   //find all the anatomy tags that have been set as reaction models
   map<int, int> ridxFromCellType;
   for (int ireaction=0; ireaction<numReactions; ++ireaction)
   {
      vector<int> anatomyCellTypes;
      objectGet(objects[ireaction], "cellTypes", anatomyCellTypes);
      for (int itag=0; itag<anatomyCellTypes.size(); ++itag)
      {
         if (ridxFromCellType.find(anatomyCellTypes[itag]) != ridxFromCellType.end())
         {
            assert(0 && "Duplicate cellTypes within the reaction models");
         }
         ridxFromCellType[anatomyCellTypes[itag]] = ireaction;
         allCellTypes_.insert(anatomyCellTypes[itag]);
      }
   }

   
   //go through the anatomy looking for cell types with no reaction
   {
      set<int> missingReactions;
      for (int icell=0; icell<cellTypes.size(); icell++)
      {
         if (cellTypes[icell] == -1) { continue; }
         if (allCellTypes_.find(cellTypes[icell]) == allCellTypes_.end()) {
            missingReactions.insert(cellTypes[icell]);
         }
      }
      for (set<int>::iterator iter=missingReactions.begin();
           iter != missingReactions.end();
           ++iter)
      {
         printf("Found cell type %d, which isn't claimed by any reaction model.\n", *iter);
      }
      assert(missingReactions.empty());
   }
      
   //count how many of each we have.
   vector<int> countFromRidx(numReactions);
   for (int ireaction=0; ireaction<numReactions; ++ireaction)
   {
      countFromRidx[ireaction] = 0;
   }
   for (int icell=0; icell<cellTypes.size(); ++icell)
   {
      if (ridxFromCellType.find(cellTypes[icell]) != ridxFromCellType.end())
      {
         countFromRidx[ridxFromCellType[cellTypes[icell]]]++;
      }
   }

   //Go through and create the extents
   extents_.resize(numReactions+1);
   extents_[0] = 0;
   for (int ireaction=0; ireaction<numReactions; ++ireaction)
   {
      int localSize = countFromRidx[ireaction];
      extents_[ireaction+1] = extents_[ireaction]+localSize;
   }

   //build up the index
   {
      EindexFromIindex_.resize(extents_[numReactions]);
      IindexFromEindex_.resize(cellTypes.size());
      wo_array_ptr<int> EindexFromIindexArray = EindexFromIindex_.useOn(CPU);
      
      vector<int> reactionCursor = extents_;
      for (int icell=0; icell<cellTypes.size(); ++icell)
      {
         int Eindex = icell;
         int Iindex = -1;
         auto ridx_iter = ridxFromCellType.find(cellTypes[icell]);
         if (ridx_iter != ridxFromCellType.end())
         {
            Iindex = reactionCursor[ridx_iter->second]++;
            EindexFromIindexArray[Iindex] = icell;
         }
         IindexFromEindex_[Eindex] = Iindex;
      }
      //finish off the index
      for (int ireaction=0; ireaction<numReactions; ++ireaction)
      {
         assert(reactionCursor[ireaction] == extents_[ireaction+1]);
      }
   }
   
   //create the reaction objects
   reactions_.resize(numReactions);
   for (int ireaction=0; ireaction<numReactions; ++ireaction)
   {
      int localSize = countFromRidx[ireaction];
      reactions_[ireaction] = reactionFactory(objectNameFromRidx_[ireaction], dt, localSize, group);
   }

   //Ok, now we've created the reaction objects.  Now we need to
   //figure out the state variables will map.

   //First, collect all the info about the state variables.
   vector<vector<string> > subUnitsFromType(numTypes);
   vector<vector<string> > subVarnamesFromType(numTypes);
   vector<vector<int> > subHandlesFromType(numTypes);
   for (int itype=0; itype<numTypes; ++itype)
   {
      //find a reaction object for this type.
      Reaction* thisReaction=NULL;
      for (int ireaction=0; ireaction<numReactions; ++ireaction)
      {
         if (typeFromRidx_[ireaction] == itype)
         {
            thisReaction = reactions_[ireaction];
            break;
         }
      }
      assert(thisReaction != NULL);
      
      //query the state variable information.
      thisReaction->getCheckpointInfo(subVarnamesFromType[itype], subUnitsFromType[itype]);
      subHandlesFromType[itype] = thisReaction->getVarHandle(subVarnamesFromType[itype]);
   }
   
   //find all the subVarnames with more than one baseType
   set<string> subVarnamesWithMultipleBases;
   {
      map<string, string> baseUnitFromSubVarname;
      for (int itype=0; itype<numTypes; ++itype)
      {
         const vector<string>& subVarnames(subVarnamesFromType[itype]);
         for (int isubVarname=0; isubVarname<subVarnames.size(); ++isubVarname)
         {
            const string& thisSubVarname(subVarnames[isubVarname]);
            const string& thisSubUnit(subUnitsFromType[itype][isubVarname]);
            //does this name has a subtype assigned to it?
            if (baseUnitFromSubVarname.find(thisSubVarname) == baseUnitFromSubVarname.end())
            {
               //if a type hasn't been identified yet for this guy, add it here.
               baseUnitFromSubVarname[thisSubVarname] = thisSubUnit;
            }
            else
            {
               //otherwise, if the units don't match this guy has multiple bases.
               if (units_check(baseUnitFromSubVarname[thisSubVarname].c_str(), thisSubUnit.c_str()))
               {
                  subVarnamesWithMultipleBases.insert(thisSubVarname);
               }
            }
         }
      }
   }
   
   //set up the rename structures
   //if a subvar has multiple bases, rename it everywhere.  Otherwise, pass it through.
   //collect up all the unique names and the bases.
   //sets unitFromHandle_
   //sets handleFromVarname_
   //sets subHandleInfoFromTypeAndHandle_
   subHandleInfoFromTypeAndHandle_.resize(numTypes);
   for (int itype=0; itype<numTypes; ++itype)
   {
      for (int isubVarname=0; isubVarname<subVarnamesFromType[itype].size(); ++isubVarname)
      {
         const string& thisSubVarname(subVarnamesFromType[itype][isubVarname]);
         const string& thisSubUnit(subUnitsFromType[itype][isubVarname]);
         const int& thisSubHandle(subHandlesFromType[itype][isubVarname]);
         string newName = thisSubVarname;
         if (subVarnamesWithMultipleBases.find(thisSubVarname) != subVarnamesWithMultipleBases.end())
         {
            newName = methodNameFromType_[itype]+"."+newName;
         }

         int thisHandle;
         if (handleFromVarname_.find(newName) == handleFromVarname_.end())
         {
            thisHandle = handleFromVarname_.size();
            handleFromVarname_[newName] = thisHandle;
            unitFromHandle_.push_back(thisSubUnit);
         }
         else
         {
            thisHandle = handleFromVarname_[newName];
         }
         assert(thisHandle == handleFromVarname_[newName]);
         assert(handleFromVarname_.size() == unitFromHandle_.size());
         assert(units_check(thisSubUnit.c_str(),unitFromHandle_[thisHandle].c_str()) == 0);
         
         subHandleInfoFromTypeAndHandle_[itype][thisHandle].first = thisSubHandle;
         subHandleInfoFromTypeAndHandle_[itype][thisHandle].second = units_convert(1,unitFromHandle_[thisHandle].c_str(),thisSubUnit.c_str());
      }
   }
}

std::string ReactionManager::stateDescription() const {
   assert(methodNameFromType_.size() >= 1);
   string retval = methodNameFromType_[0];
   for (int itype=1; itype<methodNameFromType_.size(); ++itype) {
      retval += "+";
      retval += methodNameFromType_[itype];
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
   if (handleFromVarname_.find(varName) == handleFromVarname_.end())
   {
      return -1;
   }
   else
   {
      return handleFromVarname_.find(varName)->second;
   }
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
      int subCell = IindexFromEindex_[iCell]-extents_[ridx];
      reactions_[ridx]->setValue(subCell, subHandle, value/myUnitFromTheirUnit);
   }
}
double ReactionManager::getValue(int iCell, int varHandle) const
{
   int ridx = getRidxFromCell(iCell);
   int subHandle;
   double myUnitFromTheirUnit;
   if (subUsesHandle(ridx, varHandle, subHandle, myUnitFromTheirUnit)) {
      int subCell = IindexFromEindex_[iCell]-extents_[ridx];
      return myUnitFromTheirUnit*reactions_[ridx]->getValue(subCell, subHandle);
   } else {
      return numeric_limits<double>::quiet_NaN();
   }
}
void ReactionManager::getValue(int iCell,
                               ro_array_ptr<int> handle,
                               wo_array_ptr<double> value) const
{
   for (unsigned ii=0; ii<handle.size(); ++ii)
      value[ii] = getValue(iCell, handle[ii]);
}
const std::string ReactionManager::getUnit(const std::string& varName) const
{
   return unitFromHandle_[getVarHandle(varName)];
}

vector<int> ReactionManager::allCellTypes() const
{
   vector<int> retVal(allCellTypes_.size());
   copy(allCellTypes_.begin(), allCellTypes_.end(), retVal.begin());
   return retVal;
}
   
int ReactionManager::getRidxFromCell(const int Eindex) const {
   int Iindex = IindexFromEindex_[Eindex];
   assert(Iindex != -1);
   //binary search to find the proper Ridx.
   int begin=0;
   int end=extents_.size();
   while (begin+1 < end)
   {
      int mid=(begin+end)/2;
      if (extents_[mid] == Iindex)
      {
         return mid;
      }
      else if (extents_[mid] < Iindex)
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
   const map<int, pair<int, double> >& subHandleInfoFromHandle(subHandleInfoFromTypeAndHandle_[typeFromRidx_[ridx]]);
   
   if (subHandleInfoFromHandle.find(handle) == subHandleInfoFromHandle.end())
   {
      return false;
   }

   subHandle = subHandleInfoFromHandle.find(handle)->second.first;
   myUnitFromTheirUnit = subHandleInfoFromHandle.find(handle)->second.second;
   return true;
}

