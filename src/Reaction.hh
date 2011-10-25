#ifndef REACTION_HH
#define REACTION_HH

#include <vector>

class Reaction
{
 public:
   virtual void
   calc(double dt, std::vector<double>& Vm, std::vector<double>& iStim) = 0;
   
};

#endif
