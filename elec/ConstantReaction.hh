#ifndef CONSTANTREACTION_HH
#define CONSTANTREACTION_HH

#include "Reaction.hh"
#include "ConstantModel.hh"


class ConstantReaction : public Reaction
{
 public:
   ConstantReaction(const Anatomy& anatomy,
                 const std::vector<double>& eta,
                 const SymmetricTensor& sigma1,
                 const SymmetricTensor& sigma2,
                 const SymmetricTensor& sigma3,
                 const double alpha,
                 const double beta,
                 const double gamma,
                 const int printRate);
   ~ConstantReaction(){};
   std::string methodName() const {return "Constant";}

   virtual void calc(double dt,
                     ro_larray_ptr<double> Vm,
                     ro_larray_ptr<double> iStim,
                     wo_larray_ptr<double> dVm);
   virtual void initializeMembraneVoltage(wo_larray_ptr<double> Vm);
   double getValue(int iCell, int handle) const;
   int getVarHandle(const std::string& varName) const;

 private:
   const Anatomy& anatomy_;
   ConstantModel* cellModel_;
   int printRate_;
};

#endif
