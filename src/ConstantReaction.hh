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

   void calc(double dt,
             const VectorDouble32& Vm,
             const std::vector<double>& iStim,
             VectorDouble32& dVm);
   void initializeMembraneVoltage(VectorDouble32& Vm);
   double getValue(int iCell, int handle) const;
   int getVarHandle(const string& varName) const;

 private:
   const Anatomy& anatomy_;
   ConstantModel* cellModel_;
   int printRate_;

   void compareWithExactSol(const VectorDouble32& Vm)const;
};

#endif
