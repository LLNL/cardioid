#include "Reaction.hh"
#include "object.h"
#include "TransportCoordinator.hh"
#include <vector>
#include <string>

class Anatomy;

//#define SIMD_WIDTH 2

namespace SimpleTT06 
{

   
struct State 
{
   //EDIT_STATE
   double *data;
   
   // mov nCells to the end
   int nCells;
   int nStates;
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
   TransportCoordinator<State> stateTransport_;
};

}

namespace scanReaction 
{
   Reaction* scanSimpleTT06(OBJECT* obj, const Anatomy& anatomy);
}
