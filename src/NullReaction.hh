#ifndef NULL_REACTION_HH
#define NULL_REACTION_HH

#include "Reaction.hh"

class NullReaction : public Reaction
{
 public:
   
   NullReaction(){};
   std::string methodName() const {return "null";}

   void calc(double dt,
             const std::vector<double>& Vm,
             const std::vector<double>& iStim,
             std::vector<double>& dVm){};
   void initializeMembraneVoltage(std::vector<double>& Vm){};

 private:

};

#endif
