#ifndef TT04_CELLML_REACTION_HH
#define TT04_CELLML_REACTION_HH

#include "Reaction.hh"
class Anatomy;
class TT04_CellML;
class TT04_CellMLState;

class TT04_CellML_Reaction : public Reaction
{
 public:
   enum IntegratorType {forwardEuler, rushLarsen};
   
   TT04_CellML_Reaction(const int numPoints, const int ttType, IntegratorType integrator);
   std::string methodName() const {return "TT04_CellML";}
   // copy constructor and assignment operator intentionally
   // left unimplemented.
   TT04_CellML_Reaction(const TT04_CellML_Reaction&);
   TT04_CellML_Reaction& operator=(const TT04_CellML_Reaction&);
   ~TT04_CellML_Reaction();

   void calc(double dt,
             const VectorDouble32& Vm,
             const std::vector<double>& iStim,
             VectorDouble32& dVm);
   void initializeMembraneVoltage(VectorDouble32& Vm);

 private:

   void forwardEulerIntegrator(double dt, const VectorDouble32& Vm,
      const std::vector<double>& iStim, VectorDouble32& dVm);
   void rushLarsenIntegrator(double dt, const VectorDouble32& Vm,
      const std::vector<double>& iStim, VectorDouble32& dVm);

   unsigned nCells_;
   IntegratorType integrator_;
   
   std::vector<TT04_CellML*>     cellModel_;
   std::vector<TT04_CellMLState> s_;
};

#endif
