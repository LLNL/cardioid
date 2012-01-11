#ifndef TT04Dev_REACTION_HH
#define TT04Dev_REACTION_HH

#include "Reaction.hh"
class Anatomy;
class TT04Dev;


class TT04Dev_Reaction : public Reaction
{
 public:
   TT04Dev_Reaction(const Anatomy& anatomy);
   
   void calc(double dt,
             const std::vector<double>& Vm,
             const std::vector<double>& iStim,
             std::vector<double>& dVm);
   void initializeMembraneVoltage(std::vector<double>& Vm);

 private:
   unsigned nCells_;
   std::vector<int>     ttType_; // maps cellType to ttType
   std::vector<TT04Dev> cellModel_;

};

#endif
