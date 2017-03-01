#include "Reaction.hh"
#include "object.h"
#include <vector>
class Anatomy;

namespace SimpleTT06 
{

struct State 
{
   //EDIT_STATE
   double f2Gate;
   double fGate;
   double dGate;
   double mGate;
   double jGate;
   double hGate;
   double rGate;
   double sGate;
   double Xr1Gate;
   double Xr2Gate;
   double XsGate;
   double jLGate;
   double Na_i;
   double Ca_i;
   double Ca_ss;
   double Ca_sr;
   double fCass;
   double dVK_i;
   double R_prime;
};

struct PerCellFlags 
{
   //EDIT_PERCELL_FLAGS
};

struct PerCellParameters
{
   //EDIT_PERCELL_PARAMETERS
};

class ThisReaction : public Reaction
{
 public:
   ThisReaction(const Anatomy& anatomy);
   std::string methodName() const;
   
   void calc(double dt,
             const VectorDouble32& Vm,
             const std::vector<double>& iStim,
             VectorDouble32& dVm);
   void initializeMembraneVoltage(VectorDouble32& Vm);
   virtual void getCheckpointInfo(std::vector<std::string>& fieldNames,
                                  std::vector<std::string>& fieldUnits) const;
   virtual int getVarHandle(const std::string& varName) const;
   virtual void setValue(int iCell, int varHandle, double value);
   virtual double getValue(int iCell, int varHandle) const;
   virtual const std::string getUnit(const std::string& varName) const;

 public:
   //constant flags
   //EDIT_FLAGS

   //EDIT_PARAMETERS

   //per-cell flags
   std::vector<PerCellFlags> perCellFlags_;
   std::vector<PerCellParameters> perCellParameters_;

 private:
   unsigned nCells_;
   std::vector<State> state_;
   bool isCPUValid_;
   bool isDeviceValid_;
   void updateDevice() const;
   void updateHost() const;
};

}

namespace scanReaction 
{
   Reaction* scanSimpleTT06(OBJECT* obj, const Anatomy& anatomy);
}
