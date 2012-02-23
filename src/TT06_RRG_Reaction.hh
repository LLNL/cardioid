#ifndef TT06_RRG_REACTION_HH
#define TT06_RRG_REACTION_HH

#include "Reaction.hh"
class Anatomy;
class TT06_RRG;
class BucketOfBits;

class TT06_RRG_Reaction : public Reaction
{
 public:
   
   TT06_RRG_Reaction(const Anatomy& anatomy);
   std::string methodName() const {return "TT06_RRG";}

   void calc(double dt,
             const std::vector<double>& Vm,
             const std::vector<double>& iStim,
             std::vector<double>& dVm);
   void initializeMembraneVoltage(std::vector<double>& Vm);

   void loadState(const BucketOfBits& data);

   void getCheckpointInfo(std::vector<std::string>& fieldNames,
                          std::vector<std::string>& fieldUnits) const;
   std::vector<int> getHandle(const std::vector<std::string>& varName) const;
   void getValue(int iCell,
                 const std::vector<int>& handle,
                 std::vector<double>& value) const;
   
   
 private:

   std::vector<int>      ttType_; // maps cellType to ttType
   std::vector<TT06_RRG> cells_;
};

#endif
