#ifndef NULL_REACTION_HH
#define NULL_REACTION_HH

#include "Reaction.hh"

struct NullReactionParms
{
   double initialVoltage;
};


class NullReaction : public Reaction
{
 public:
   
   NullReaction(const NullReactionParms& parms)
   : V0_(parms.initialVoltage)
   {};
   std::string methodName() const {return "null";}

   void calc(double dt,
             const VectorDouble32& Vm,
             const std::vector<double>& iStim,
             VectorDouble32& dVm){};
   void initializeMembraneVoltage(VectorDouble32& Vm)
   {
      for (unsigned ii=0; ii<Vm.size(); ++ii)
         Vm[ii] = V0_;
   };

 private:
   double V0_;

};

#endif
