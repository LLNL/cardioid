#ifndef TT04_CELLML_REACTION_HH
#define TT04_CELLML_REACTION_HH

#include "Reaction.hh"
class Anatomy;
class TT04_CellML;
class TT04_CellMLState;

class TT04_CellML_Reaction : public Reaction
{
 public:
   enum IntegratorType {forwardEuler, rushLarson};
   
   TT04_CellML_Reaction(const Anatomy& anatomy, IntegratorType integrator);
   // copy constructor and assignment operator intentionally
   // left unimplemented.
   TT04_CellML_Reaction(const TT04_CellML_Reaction&);
   TT04_CellML_Reaction& operator=(const TT04_CellML_Reaction&);
   ~TT04_CellML_Reaction();

   void calc(double dt,
             const std::vector<double>& Vm,
             const std::vector<double>& iStim,
             std::vector<double>& dVm);

 private:

   void forwardEulerIntegrator(double dt, const std::vector<double>& Vm,
      const std::vector<double>& iStim, std::vector<double>& dVm);
   void rushLarsonIntegrator(double dt, const std::vector<double>& Vm,
      const std::vector<double>& iStim, std::vector<double>& dVm);

   unsigned nCells_;
   IntegratorType integrator_;
   
   std::vector<int>              ttType_; // maps cellType to ttType
   std::vector<TT04_CellML*>     cellModel_;
   std::vector<TT04_CellMLState> s_;
};

#endif
