/**

   How to convert this code to work for any other model:

   - Search/Replace the model name with your own specific string in the header and source files
   - Add your own code to EDIT_FLAGS and EDIT_PARAMETERS
   - Add your own code to EDIT_PERCELL_FLAGS and EDIT_PERCELL_PARAMETERS
   - Add your own states to EDIT_STATE
   - Add your computation code to the main calc routine, copy pasting frmo matlab.
   
 */


#include "SimpleTT06.hh"
#include "object_cc.hh"
#include "Anatomy.hh"
//#include "jitifyDecl.hpp"
#include <cmath>

#define pi 3.141592653589793238462643383279502884197169399375105820974944592307816406286

#define fCassForm TT06

using namespace std;

namespace scanReaction 
{

#define setDefault(name, value) objectGet(obj, #name, reaction->name, #value)
   
   Reaction* scanSimpleTT06(OBJECT* obj, const int numPoints)
   {
      SimpleTT06::ThisReaction* reaction = new SimpleTT06::ThisReaction(numPoints);

      //override the defaults
      //EDIT_FLAGS

      //EDIT_PARAMETERS
      
      return reaction;
   }
#undef setDefault

}

namespace SimpleTT06 
{
   
string ThisReaction::methodName() const
{
   return "SimpleTT06";
}
const char* varNames[] = 
{
   //EDIT_STATE
   "f2Gate",
   "fGate",
   "dGate",
   "mGate",
   "jGate",
   "hGate",
   "rGate",
   "sGate",
   "Xr1Gate",
   "Xr2Gate",
   "XsGate",
   "jLGate",
   "Na_i",
   "Ca_i",
   "Ca_ss",
   "Ca_sr",
   "fCass",
   "dVK_i",
   "R_prime"
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
   return -1;
}

void assertStateOrderAndVarNamesAgree(void)
{
   //State s;
//#define checkVarOrder(x) assert(reinterpret_cast<double*>(&s)+getVarOffset(#x) == &s . x)

   //int STATIC_ASSERT_checkAllDouble[(NUMVARS == sizeof(s)/sizeof(double))? 1: 0];

   //EDIT_STATE
/*   checkVarOrder(f2Gate);
   checkVarOrder(fGate);
   checkVarOrder(dGate);
   checkVarOrder(mGate);
   checkVarOrder(jGate);
   checkVarOrder(hGate);
   checkVarOrder(rGate);
   checkVarOrder(sGate);
   checkVarOrder(Xr1Gate);
   checkVarOrder(Xr2Gate);
   checkVarOrder(XsGate);
   checkVarOrder(jLGate);
   checkVarOrder(Na_i);
   checkVarOrder(Ca_i);
   checkVarOrder(Ca_ss);
   checkVarOrder(Ca_sr);
   checkVarOrder(fCass);
   checkVarOrder(dVK_i);
   checkVarOrder(R_prime);
*/
   }

enum StateOffset {
   f2Gate_off,
   fGate_off,
   dGate_off,
   mGate_off,
   jGate_off,
   hGate_off,
   rGate_off,
   sGate_off,
   Xr1Gate_off,
   Xr2Gate_off,
   XsGate_off,
   jLGate_off,
   Na_i_off,
   Ca_i_off,
   Ca_ss_off,
   Ca_sr_off,
   fCass_off,
   dVK_i_off,
   R_prime_off,
   NUMSTATES
};
   
ThisReaction::ThisReaction(const int numPoints)
: nCells_(numPoints)
{   
   stateTransport_.setup(PinnedVector<double>(nCells_*NUMSTATES));
   
   assertStateOrderAndVarNamesAgree();
   perCellFlags_.resize(nCells_);
   perCellParameters_.resize(nCells_);   
}
   
   void actualCalcJitify(const double dt, OnDevice<ConstArrayView<double>> Vm, OnDevice<ConstArrayView<double>> iStim, OnDevice<ArrayView<double>> dVm, OnDevice<ArrayView<double>> stateData);
   
void ThisReaction::calc(double dt,
                const Managed<ArrayView<double>> Vm,
                const Managed<ArrayView<double>> iStim,
                Managed<ArrayView<double>> dVm)
{   
   actualCalcJitify(dt,
                    Vm,
                    iStim,
                    dVm,
                    stateTransport_);
}
   
void ThisReaction::initializeMembraneVoltage(ArrayView<double> Vm)
{
   assert(Vm.size() >= nCells_);
   const double initVm = -86.709;
   Vm.assign(Vm.size(), initVm);
   //State initState;
   //EDIT_STATE

   ArrayView<double> stateData = stateTransport_;

   for (int ii=0; ii<nCells_; ii++)
   {
      const double pcnst_2 = 96485.3415;
      const double pcnst_3 = 0.185;
      const double pcnst_4 = 0.016404;
      const double c9 = -pcnst_3/(pcnst_4*pcnst_2);
      const double K_i = 138.4;
      stateData[dVK_i_off*nCells_+ii] = K_i/c9+initVm;
      stateData[Na_i_off*nCells_+ii]     =10.355;
      stateData[Ca_i_off*nCells_+ii]     =0.00013;
      stateData[Ca_ss_off*nCells_+ii]    =0.00036 ;
      stateData[Ca_sr_off*nCells_+ii]    =3.715   ;
      stateData[R_prime_off*nCells_+ii]  =0.9068  ;
      stateData[fCass_off*nCells_+ii]    =0.9953  ;
      stateData[Xr1Gate_off*nCells_+ii]  =0.00448 ;
      stateData[Xr2Gate_off*nCells_+ii]  =0.476   ;
      stateData[XsGate_off*nCells_+ii]   =0.0087  ;
      stateData[mGate_off*nCells_+ii]    =0.00155 ;
      stateData[hGate_off*nCells_+ii]    =0.7573  ;
      stateData[jGate_off*nCells_+ii]    =0.7225  ;
      stateData[rGate_off*nCells_+ii]    =2.235e-8;
      stateData[dGate_off*nCells_+ii]    =3.164e-5;
      stateData[fGate_off*nCells_+ii]    =0.8009  ;
      stateData[f2Gate_off*nCells_+ii]   =0.9778  ;
      stateData[sGate_off*nCells_+ii]    =0.3212  ;
      stateData[jLGate_off*nCells_+ii]   =0.066   ;
   }
   
   //state.assign(state.size(), initState);
   
}

const string ThisReaction::getUnit(const std::string& varName) const
{
   //deliberatly broken for now, if this code still is being used past 2016-11-01 something has gone wrong.
   return "1";
}

#define HANDLE_OFFSET 1000
int ThisReaction::getVarHandle(const std::string& varName) const
{
   int retVal = getVarOffset(varName);
   if (retVal >= 0)
   {
      retVal += HANDLE_OFFSET;
   }
   return retVal;
}

void ThisReaction::setValue(int iCell, int varHandle, double value) 
{
   ArrayView<double> stateData = stateTransport_;
   int var=varHandle-HANDLE_OFFSET;

   stateData[var*nCells_+iCell] = value;
}


double ThisReaction::getValue(int iCell, int varHandle) const
{
   ConstArrayView<double> stateData = stateTransport_;
   int var=varHandle-HANDLE_OFFSET;
  
   return stateData[var*nCells_+iCell]; 
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
