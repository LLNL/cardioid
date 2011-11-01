#ifndef REACTION_HH
#define REACTION_HH

#include <vector>

class Reaction
{
 public:
   virtual void calc(double dt,
		     const std::vector<double>& Vm,
		     const std::vector<double>& iStim,
		     std::vector<double>& dVm) = 0;
   
};

#endif
