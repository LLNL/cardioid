#ifndef OHaraRudy_REACTION_HH
#define OHaraRudy_REACTION_HH

#include "Reaction.hh"
class Anatomy;
class OHaraRudy;
class BucketOfBits;

class OHaraRudy_Reaction : public Reaction
{
 public:
   
   OHaraRudy_Reaction(const Anatomy& anatomy);
   std::string methodName() const {return "OHaraRudy";}

   void calc(double dt,
             const VectorDouble32& Vm,
             const std::vector<double>& iStim,
             VectorDouble32& dVm);
   void initializeMembraneVoltage(VectorDouble32& Vm);

   /** Functions needed for checkpoint/restart */
   void getCheckpointInfo(std::vector<std::string>& fieldNames,
                          std::vector<std::string>& fieldUnits) const;
   int getVarHandle(const std::string& varName) const;
   void setValue(int iCell, int varHandle, double value);
   double getValue(int iCell, int varHandle) const;
   void getValue(int iCell,
                 const std::vector<int>& handle,
                 std::vector<double>& value) const;
   const std::string getUnit(const std::string& varName) const;
   
   
 private:

   int indexS_;
   std::vector<int>      ttType_; // maps cellType to ttType
   std::vector<OHaraRudy> cells_;
};

#endif
