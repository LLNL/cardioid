/**

   How to convert this code to work for any other model:

   - Search/Replace the model name with your own specific string in the header and source files
   - Add your own code to EDIT_FLAGS and EDIT_PARAMETERS
   - Add your own code to EDIT_PERCELL_FLAGS and EDIT_PERCELL_PARAMETERS
   - Add your own states to EDIT_STATE
   - Add your computation code to the main calc routine, copy pasting frmo matlab.
   
 */


#include "Passive.hh"
#include "object_cc.hh"
#include "Anatomy.hh"
#include <cmath>

#define pi 3.141592653589793238462643383279502884197169399375105820974944592307816406286

using namespace std;

namespace scanReaction 
{

#define setDefault(name, value) objectGet(obj, #name, reaction->name, #value)
   
   Reaction* scanPassive(OBJECT* obj, const int numPoints)
   {
      Passive::ThisReaction* reaction = new Passive::ThisReaction(numPoints);

      //override the defaults
      //EDIT_FLAGS

      //EDIT_PARAMETERS
      setDefault(G, 6.0643e-4);    // [mS/uF] 3
      setDefault(E_R, -85.);
      
      return reaction;
   }
#undef setDefault

}

namespace Passive 
{

inline double pow(const double x, const int p)
{
   double ret=1;
   if (p > 0) 
   {
      for (int ii=0; ii<p; ii++) 
      {
         ret *= x;
      }
   }
   else
   {
      for (int ii=0; ii<-p; ii++) 
      {
         ret /= x;
      }
   }
   return ret;
}
   
string ThisReaction::methodName() const
{
   return "Passive";
}
const char* varNames[] = 
{
   //EDIT_STATE
};
#define NUMVARS (sizeof(varNames)/sizeof(char*))

int getVarOffset(const std::string& varName)
{
   for (int ivar=0; ivar<NUMVARS; ivar++) 
   {
      if (varNames[ivar] == varName) 
      {
         return ivar;
      }
   }
   assert(0 && "Control should never get here.");
   return -1;
}

void assertStateOrderAndVarNamesAgree(void)
{
   State s;
#define checkVarOrder(x) assert(reinterpret_cast<double*>(&s)+getVarOffset(#x) == &s . x)

   int STATIC_ASSERT_checkAllDouble[(NUMVARS == sizeof(s)/sizeof(double))? 1: 0];

   //EDIT_STATE
}

   
ThisReaction::ThisReaction(const int numPoints)
: nCells_(numPoints)
{
   assertStateOrderAndVarNamesAgree();
   state_.resize(nCells_);
   perCellFlags_.resize(nCells_);
   perCellParameters_.resize(nCells_);
}

void ThisReaction::calc(double dt, const VectorDouble32& Vm,
                       const vector<double>& iStim , VectorDouble32& dVm)
{
   for (unsigned ii=0; ii<nCells_; ++ii)
   {
      const double v = Vm[ii];
      dVm[ii] = G*(v-E_R);      
   }
}

void ThisReaction::initializeMembraneVoltage(VectorDouble32& Vm)
{
   assert(Vm.size() >= nCells_);
   Vm.assign(Vm.size(), -85.0);
   State initState;
   //EDIT_STATE
   state_.resize(nCells_);
   state_.assign(state_.size(), initState);

}

const string ThisReaction::getUnit(const std::string& varName) const
{
   //deliberatly broken for now, if this code still is being used past 2016-11-01 something has gone wrong.
   return "1";
}

#define HANDLE_OFFSET 1000
int ThisReaction::getVarHandle(const std::string& varName) const
{
   return getVarOffset(varName)+HANDLE_OFFSET;
}

void ThisReaction::setValue(int iCell, int varHandle, double value) 
{
   reinterpret_cast<double*>(&state_[iCell])[varHandle-HANDLE_OFFSET] = value;
}


double ThisReaction::getValue(int iCell, int varHandle) const
{
   return reinterpret_cast<const double*>(&state_[iCell])[varHandle-HANDLE_OFFSET];
}

void ThisReaction::getCheckpointInfo(vector<string>& fieldNames,
                                     vector<string>& fieldUnits) const
{
   fieldNames.resize(NUMVARS);
   fieldUnits.resize(NUMVARS);
   for (int ivar=0; ivar<NUMVARS; ivar++) 
   {
      fieldNames[ivar] = varNames[ivar];
      fieldUnits[ivar] = getUnit(fieldNames[ivar]);
   }
}

}
