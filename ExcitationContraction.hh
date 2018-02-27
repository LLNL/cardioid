#pragma once

#include <vector>
#include <string>

class XXX {
 public:
   virtual int getHandle(const std::string& name) const = 0;
   virtual std::string getUnit(const int handle) const = 0;

   virtual std::vector<std::string> getVarGroup(const std::string groupName) const = 0;

   virtual double get(const int varHandle) const = 0;
   virtual void   set(const int varHandle, const double value) = 0;

   virtual double get(const int varHandle, const int iCell) const = 0;
   virtual void   set(const int varHandle, const int iCell, const double value) = 0;

   virtual double get(const int varHandle, const int iCell, const double* inputs[]) const = 0;
   
   std::vector<int> getHandles(const std::vector<std::string>& names) const;
   bool getOrder(const std::string type, const std::vector<std::string>& names, std::vector<int>& permutation) const;
   virtual ~XXX() {};
};


class ExcitationContraction : public XXX {
 public:
   virtual std::string name() const = 0;

   virtual void resetDefaultParameters() = 0;
   virtual void initialize(const double *const inputs[]) = 0;
   virtual void tryTimestep(const double dt, const double *const inputs[]) = 0;
   virtual void outputForNextTimestep(const double dt, const double* const inputs[], double* const outputs[]) = 0;
   virtual void commitTimestep() = 0;

   virtual ~ExcitationContraction() {};

;
};


