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
             const VectorDouble32& Vm,
             const std::vector<double>& iStim,
             VectorDouble32& dVm);
   void initializeMembraneVoltage(VectorDouble32& Vm);

 private:
   unsigned nCells_;
   std::vector<int>     ttType_; // maps cellType to ttType
   std::vector<TT04Dev> cellModel_;

};

#endif
