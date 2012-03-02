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

   /** Functions needed for checkpoint/restart */
   void getCheckpointInfo(std::vector<std::string>& fieldNames,
                          std::vector<std::string>& fieldUnits) const;
   int getVarHandle(const std::string& varName) const;
   std::vector<int> getVarHandle(const std::vector<std::string>& varName) const;
   void setValue(int iCell, int varHandle, double value);
   double getValue(int iCell, int varHandle) const;
   void getValue(int iCell,
                 const std::vector<int>& handle,
                 std::vector<double>& value) const;
   const std::string getUnit(const std::string& varName) const;
   
   
 private:

   std::vector<int>      ttType_; // maps cellType to ttType
   std::vector<TT06_RRG> cells_;
};

#endif
