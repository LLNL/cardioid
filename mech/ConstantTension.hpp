#include "ExcitationContraction.hpp"
#include <vector>
#include <string>

namespace ConstantTension
{

class ThisModel : public ExcitationContraction
{
 public:
   ThisModel(const int numPoints);
   virtual std::string name() const;
   virtual int getHandle(const std::string& varName) const;
   virtual std::string getUnit(const int varHandle) const;

   virtual std::vector<std::string> getVarGroup(const std::string type) const;
   
   virtual double get(const int varHandle) const;
   virtual void   set(const int varHandle, const double value);

   virtual double get(const int varHandle, const int iCell) const;
   virtual void   set(const int varHandle, const int iCell, const double value);

   virtual double get(const int varHandle, const int iCell, const double* inputs[]) const;

   virtual void resetDefaultParameters();
   virtual void initialize(const double *const inputs[]);
   virtual void tryTimestep(const double dt, const double *const inputs[]);
   virtual void outputForNextTimestep(const double dt, const double* const inputs[], double* const outputs[]);
   virtual void commitTimestep();

 private:
   int _numPoints;
   double _dt;

   //PARAMETERS
   double usedTension;
};

}
