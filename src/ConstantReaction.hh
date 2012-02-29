#ifndef CONSTANTREACTION_HH
#define CONSTANTREACTION_HH

#include "Reaction.hh"
#include "ConstantModel.hh"


class ConstantReaction : public Reaction
{
 public:
   ConstantReaction(const Anatomy& anatomy,
                 const double eta[3],
                 const SymmetricTensor& sigma1,
                 const SymmetricTensor& sigma2,
                 const SymmetricTensor& sigma3,
                 const double alpha,
                 const double beta,
                 const double gamma);
   ~ConstantReaction(){};
   std::string methodName() const {return "Constant";}

   void calc(double dt,
             const std::vector<double>& Vm,
             const std::vector<double>& iStim,
             std::vector<double>& dVm);
   void initializeMembraneVoltage(std::vector<double>& Vm);

 private:
   const Anatomy& anatomy_;
   ConstantModel* cellModel_;

};

#endif
