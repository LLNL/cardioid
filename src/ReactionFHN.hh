#ifndef REACTION_FHN_HH
#define REACTION_FHN_HH

#include "Reaction.hh"
class Anatomy;

class ReactionFHN : public Reaction
{
 public:
   ReactionFHN(const Anatomy& anatomy);
   
   void calc(double dt,
             const std::vector<double>& Vm,
             const std::vector<double>& iStim,
             std::vector<double>& dVm);
   void initializeMembraneVoltage(std::vector<double>& Vm);

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
