#ifndef TT06_RRG_REACTION_HH
#define TT06_RRG_REACTION_HH

#include "Reaction.hh"
class Anatomy;
class TT06_RRG;

class TT06_RRG_Reaction : public Reaction
{
 public:
   
   TT06_RRG_Reaction(const Anatomy& anatomy);

   void calc(double dt,
             const std::vector<double>& Vm,
             const std::vector<double>& iStim,
             std::vector<double>& dVm);

 private:

   std::vector<int>      ttType_; // maps cellType to ttType
   std::vector<TT06_RRG> cells_;
};

#endif
