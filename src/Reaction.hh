#ifndef REACTION_HH
#define REACTION_HH

#include <vector>
#include <string>
#include "VectorDouble32.hh"

class BucketOfBits;

class Reaction
{
 public:
   virtual ~Reaction(){};
   virtual std::string methodName() const = 0;
   virtual void calc(double dt,
                     const VectorDouble32& Vm,
                     const std::vector<double>& iStim,
                     VectorDouble32& dVm) = 0;
   virtual void updateNonGate(double dt, const VectorDouble32& Vm, VectorDouble32& dVR) {};
   virtual void updateGate   (double dt, const VectorDouble32& Vm) {};

   /** Populates the Vm array with some sensible default initial
    * membrane voltage.  Vm will be the parallel to the local cells in
    * the anatomy that was used to create the concrete reaction class. */
   virtual void initializeMembraneVoltage(VectorDouble32& Vm) = 0;

   virtual void scaleCurrents(std::vector<double>);  

   /** Functions needed for checkpoint/restart */
   virtual void getCheckpointInfo(std::vector<std::string>& fieldNames,
                                  std::vector<std::string>& fieldUnits) const;
   virtual int getVarHandle(const std::string& varName) const;
   std::vector<int> getVarHandle(const std::vector<std::string>& varName) const;
   virtual void setValue(int iCell, int varHandle, double value);
   virtual double getValue(int iCell, int varHandle) const;
   virtual void getValue(int iCell,
                         const std::vector<int>& handle,
                         std::vector<double>& value) const;
   virtual const std::string getUnit(const std::string& varName) const;
};

#endif
