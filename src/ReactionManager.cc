
#include "ReactionManager.hh"
#include "Reaction.hh"

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
      string objectName = rxnObjectNames_[ii];
      ::initializeMembraneState(reaction, objectName, Vm);
   }
}

std::string ReactionManager::stateDescription() const {
   return "";
}


void ReactionManager::addReaction(const std::string& rxnObjectName)
{
   rxnObjectNames_.push_back(rxnObjectName);
}

void ReactionManager::create(const double dt, Anatomy& anatomy, const ThreadTeam &group, const std::vector<std::string>& scaleCurrents)
{
}

   /** Functions needed for checkpoint/restart */
void ReactionManager::getCheckpointInfo(std::vector<std::string>& fieldNames,
                                        std::vector<std::string>& fieldUnits) const
{
}
int ReactionManager::getVarHandle(const std::string& varName) const
{
   return 0;
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
}
double ReactionManager::getValue(int iCell, int varHandle) const
{
   return 0;
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
   return "1";
}
   
