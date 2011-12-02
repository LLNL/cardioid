#ifndef TT04_BB_REACTION_HH
#define TT04_BB_REACTION_HH

#include "Reaction.hh"
class Anatomy;
class IBM_tenTusscher04;


class TT04_bbReaction : public Reaction
{
 public:
   TT04_bbReaction(const Anatomy& anatomy);
   
   void calc(double dt,
             const std::vector<double>& Vm,
             const std::vector<double>& iStim,
             std::vector<double>& dVm);

 private:
   unsigned nCells_;
   std::vector<int>               ttType_; // maps cellType to ttType
   std::vector<IBM_tenTusscher04> cellModel_;

};

#endif
