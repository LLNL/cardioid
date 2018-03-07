#ifndef CONSTANTREACTION_HH
#define CONSTANTREACTION_HH

#include "Reaction.hh"
#include "ConstantModel.hh"


class ConstantReaction : public Reaction
{
 public:
   ConstantReaction(const Anatomy& anatomy,
                 const vector<double>& eta,
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
                     const Managed<ArrayView<double>> Vm,
                     const Managed<ArrayView<double>> iStim,
                     Managed<ArrayView<double>> dVm);
   virtual void initializeMembraneVoltage(ArrayView<double> Vm);
   double getValue(int iCell, int handle) const;
   int getVarHandle(const string& varName) const;

 private:
   const Anatomy& anatomy_;
   ConstantModel* cellModel_;
   int printRate_;
};

#endif
