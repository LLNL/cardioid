#ifndef TT04_REACTION_HH
#define TT04_REACTION_HH

#include "Reaction.hh"
class Anatomy;
class IBM_TT04;


class TT04Reaction : public Reaction
{
 public:
   TT04Reaction(const Anatomy& anatomy);
   
   void calc(double dt,
	     const std::vector<double>& Vm,
	     const std::vector<double>& iStim,
	     std::vector<double>& dVm);

 private:
   unsigned nCells_;
   std::vector<int>      ttType_; // maps cellType to ttType
   std::vector<IBM_TT04> cellModel_;

};

#endif
