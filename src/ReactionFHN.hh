#ifndef REACTION_FHN_HH
#define REACTION_FHN_HH

#include "Reaction.hh"
class Anatomy;

class ReactionFHN : public Reaction
{
 public:
   ReactionFHN(const int numPoints);
   std::string methodName() const {return "FHN";}
   
   void calc(double dt,
             const Managed<ArrayView<double>> Vm,
             const Managed<ArrayView<double>> iStim,
             Managed<ArrayView<double>> dVm);
   void initializeMembraneVoltage(ArrayView<double> Vm);

 private:
   unsigned nCells_;
   double fhn1_;
   double fhn2_;
   double Vrest_;
   double Vthresh_;
   double Vpeak_;
   double b_;
   double c3_;

   std::vector<double> W_;
   
};

#endif
