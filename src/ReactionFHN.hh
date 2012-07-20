#ifndef REACTION_FHN_HH
#define REACTION_FHN_HH

#include "Reaction.hh"
class Anatomy;

class ReactionFHN : public Reaction
{
 public:
   ReactionFHN(const Anatomy& anatomy);
   std::string methodName() const {return "FHN";}
   
   void calc(double dt,
             const VectorDouble32& Vm,
             const std::vector<double>& iStim,
             VectorDouble32& dVm);
   void initializeMembraneVoltage(VectorDouble32& Vm);

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
