#ifndef LUMENS
#define LUMENS

#include "ExcitationContraction.hpp"
#include <vector>
#include <string>

namespace Lumens2009
{

struct State
{
   double C;
   double Lsc;
   double _partial_C;
   double _partial_Lsc;
};

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
   //double _dt;

   //FLAGS

   double celltype;
   //PARAMETERS

   //STATE
   std::vector<State> _state[2];
   std::vector<State> _partial;
   int _oldIdx;

   inline int _old() const { return _oldIdx; }
   inline int _new() const { return 1 ^ _oldIdx; }
   inline void swapOldAndNew() { _oldIdx = _new(); }
  
};

}
   

#endif
