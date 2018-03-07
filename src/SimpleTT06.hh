#include "Reaction.hh"
#include "object.h"
#include "TransportCoordinator.hh"
#include <vector>
#include <string>

class Anatomy;

//#define SIMD_WIDTH 2

namespace SimpleTT06 
{

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
   ThisReaction(const int numPoints);
   std::string methodName() const;
   
   void calc(double dt,
             const Managed<ArrayView<double>> Vm_m,
             const Managed<ArrayView<double>> iStim_m,
             Managed<ArrayView<double>> dVm_m);
   void initializeMembraneVoltage(ArrayView<double> Vm);
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
   TransportCoordinator<PinnedVector<double>> stateTransport_;
};

}

namespace scanReaction 
{
   Reaction* scanSimpleTT06(OBJECT* obj, const Anatomy& anatomy);
}
