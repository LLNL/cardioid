#ifndef NULL_REACTION_HH
#define NULL_REACTION_HH

#include "Reaction.hh"

class NullReaction : public Reaction
{
 public:
   
   NullReaction(){};

   void calc(double dt,
             const std::vector<double>& Vm,
             const std::vector<double>& iStim,
             std::vector<double>& dVm){};

 private:

};

#endif
