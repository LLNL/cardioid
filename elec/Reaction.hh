#ifndef REACTION_HH
#define REACTION_HH

#include <vector>
#include <string>
#include "lazy_array.hh"

class BucketOfBits;

class Reaction
{
 public:
   virtual ~Reaction(){};
   virtual std::string methodName() const = 0;
   virtual void calc(double dt,
                     ro_mgarray_ptr<int> indexArray,
                     ro_mgarray_ptr<double> Vm,
                     ro_mgarray_ptr<double> iStim,
                     wo_mgarray_ptr<double> dVm) = 0;
   virtual void updateNonGate(double dt, ro_mgarray_ptr<int> indexArray, ro_mgarray_ptr<double> Vm, wo_mgarray_ptr<double> dVR) {};
   virtual void updateGate   (double dt, ro_mgarray_ptr<int> indexArray, ro_mgarray_ptr<double> Vm) {};

   /** Populates the Vm array with some sensible default initial
    * membrane voltage.  Vm will be the parallel to the local cells in
    * the anatomy that was used to create the concrete reaction class. */
   virtual void initializeMembraneVoltage(ro_mgarray_ptr<int> indexArray, wo_mgarray_ptr<double> Vm) = 0;

   /** Functions needed for checkpoint/restart */
   virtual void getCheckpointInfo(std::vector<std::string>& fieldNames,
                                  std::vector<std::string>& fieldUnits) const;
   virtual int getVarHandle(const std::string& varName) const;
   std::vector<int> getVarHandle(const std::vector<std::string>& varName) const;
   virtual void setValue(int iCell, int varHandle, double value);
   virtual double getValue(int iCell, int varHandle) const;
   virtual double getValue(int iCell, int varHandle, double V) const;
   virtual void getValue(int iCell,
                         const std::vector<int>& handle,
                         std::vector<double>& value) const;
   virtual const std::string getUnit(const std::string& varName) const;
};

//! Call this instead of initializeMembraneVoltage directly.
void initializeMembraneState(Reaction* reaction, const std::string& objectName, ro_mgarray_ptr<int> indexArray, wo_mgarray_ptr<double> Vm);

#endif
