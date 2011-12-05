#ifndef TT06_CELLML_REACTION_HH
#define TT06_CELLML_REACTION_HH

#include "Reaction.hh"
class Anatomy;
class TT06_CellML;
class TT06_CellMLState;

class TT06_CellML_Reaction : public Reaction
{
 public:
   enum IntegratorType {forwardEuler, rushLarsen};
   
   TT06_CellML_Reaction(const Anatomy& anatomy, IntegratorType integrator);
   // copy constructor and assignment operator intentionally
   // left unimplemented.
   TT06_CellML_Reaction(const TT06_CellML_Reaction&);
   TT06_CellML_Reaction& operator=(const TT06_CellML_Reaction&);
   ~TT06_CellML_Reaction();

   void calc(double dt,
             const std::vector<double>& Vm,
             const std::vector<double>& iStim,
             std::vector<double>& dVm);

 private:

   void forwardEulerIntegrator(double dt, const std::vector<double>& Vm,
      const std::vector<double>& iStim, std::vector<double>& dVm);
   void rushLarsenIntegrator(double dt, const std::vector<double>& Vm,
      const std::vector<double>& iStim, std::vector<double>& dVm);

   unsigned nCells_;
   IntegratorType integrator_;
   
   std::vector<int>              ttType_; // maps cellType to ttType
   std::vector<TT06_CellML*>     cellModel_;
   std::vector<TT06_CellMLState> s_;
};

#endif
