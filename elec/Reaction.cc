#include "Reaction.hh"
#include "object_cc.hh"
#include <cassert>

using namespace std;

void initializeMembraneState(Reaction* reaction, const string& objectName, ro_mgarray_ptr<int> indexArray, wo_mgarray_ptr<double> _Vm)
{
   reaction->initializeMembraneVoltage(indexArray, _Vm);
   OBJECT* reactionObj = objectFind(objectName, "REACTION");
   assert(reactionObj != NULL);
   if (object_testforkeyword(reactionObj, "initialState"))
   {
      //Get the state object.
      string stateName;
      objectGet(reactionObj, "initialState", stateName, "");
      assert(!stateName.empty());
      OBJECT* stateObj = objectFind(stateName, "SINGLECELL");
      assert(stateObj != NULL);
      string stateMethodName;
      objectGet(stateObj, "method", stateMethodName, "");
      assert(stateMethodName == reaction->methodName());

      //read in the voltage
      if (object_testforkeyword(stateObj, "Vm"))
      {
         wo_array_ptr<double> Vm = _Vm.useOn(CPU);
         double newVoltage;
         objectGet(stateObj, "Vm", newVoltage, "", "mV");
         //set the voltage
         for (int ii=0; ii<Vm.size(); ++ii)
         {
            Vm[ii] = newVoltage;
         }
      }

      //read in the rest of the state variables
      vector<string> reactionStateNames;
      vector<string> reactionStateUnits;
      reaction->getCheckpointInfo(reactionStateNames, reactionStateUnits);
      for (int istate=0; istate<reactionStateNames.size(); istate++)
      {
         string& thisState = reactionStateNames[istate];
         string& thisUnit = reactionStateUnits[istate];
         if (object_testforkeyword(stateObj, thisState.c_str()))
         {
            double newValue;
            objectGet(stateObj, thisState, newValue, "", thisUnit);
            int handle = reaction->getVarHandle(thisState);
            for (int ii=0; ii<_Vm.size(); ++ii)
            {
               reaction->setValue(ii, handle, newValue);
            }
         }
      }
   }
}

/** Gets the names and units of all fields that should be written in a
 *  checkpoint file. */
void Reaction::getCheckpointInfo(vector<string>& fieldNames,
                                 vector<string>& fieldUnits) const
{
   // Concrete classes must provide override to enable
   // checkpoint/restart.
   assert(false);
}

/** getVarHandle must return -1 for any name that is not defined in
 *  the concrete reaction class. */
int Reaction::getVarHandle(const string& varName) const
{
   assert(false);
   return 0;
}

/** getVarHandle must return -1 for any name that is not defined in
 *  the concrete reaction class. */
vector<int> Reaction::getVarHandle(const vector<string>& varName) const
{
   vector<int> handle;
   for (unsigned ii=0; ii<varName.size(); ++ii)
      handle.push_back(getVarHandle(varName[ii]));

   return handle;
}

void Reaction::setValue(int iCell, int varHandle, double value)
{
   assert(false);
}


double Reaction::getValue(int iCell, int handle, double V) const
{
   return getValue(iCell, handle);
}

double Reaction::getValue(int iCell, int handle) const
{
   assert(false);
   return 0.0;
}

void Reaction::getValue(int iCell,
                        const vector<int>& handle,
                        vector<double>& value) const
{
   for (unsigned ii=0; ii<handle.size(); ++ii)
      value[ii] = getValue(iCell, handle[ii]);
}

const string Reaction::getUnit(const string& varName) const
{
   assert(false);
   return string();
}
