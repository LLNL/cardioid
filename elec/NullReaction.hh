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
             const Managed<ArrayView<double>> Vm,
             const Managed<ArrayView<double>> iStim,
             Managed<ArrayView<double>> dVm) {}
   void initializeMembraneVoltage(ArrayView<double> Vm)
   {
      for (unsigned ii=0; ii<Vm.size(); ++ii)
         Vm[ii] = V0_;
   };
   void getCheckpointInfo(std::vector<std::string>& fieldNames,
                          std::vector<std::string>& fieldUnits) const {
      fieldNames.clear();
      fieldUnits.clear();
   }
   virtual int getVarHandle(const std::string& varName) const {
      return -1;
   }
   
 private:
   double V0_;

};

#endif
