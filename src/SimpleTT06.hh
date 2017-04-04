#include "Reaction.hh"
#include "object.h"
#include "TransportCoordinator.hh"
#include <vector>

class Anatomy;

#define SIMD_WIDTH 2

namespace SimpleTT06 
{

   
struct alignas(sizeof(double)*SIMD_WIDTH) State 
{
   //EDIT_STATE
   double f2Gate[SIMD_WIDTH];
   double fGate[SIMD_WIDTH];
   double dGate[SIMD_WIDTH];
   double mGate[SIMD_WIDTH];
   double jGate[SIMD_WIDTH];
   double hGate[SIMD_WIDTH];
   double rGate[SIMD_WIDTH];
   double sGate[SIMD_WIDTH];
   double Xr1Gate[SIMD_WIDTH];
   double Xr2Gate[SIMD_WIDTH];
   double XsGate[SIMD_WIDTH];
   double jLGate[SIMD_WIDTH];
   double Na_i[SIMD_WIDTH];
   double Ca_i[SIMD_WIDTH];
   double Ca_ss[SIMD_WIDTH];
   double Ca_sr[SIMD_WIDTH];
   double fCass[SIMD_WIDTH];
   double dVK_i[SIMD_WIDTH];
   double R_prime[SIMD_WIDTH];
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
   TransportCoordinator<std::vector<State> > stateTransport_;
};

}

namespace scanReaction 
{
   Reaction* scanSimpleTT06(OBJECT* obj, const Anatomy& anatomy);
}
