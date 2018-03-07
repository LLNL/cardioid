#ifndef TT06_RRG_REACTION_HH
#define TT06_RRG_REACTION_HH

#include "Reaction.hh"
class Anatomy;
class TT06_RRG;
class BucketOfBits;

struct TT06_RRG_ReactionParms
{
   double Ko;
};


class TT06_RRG_Reaction : public Reaction
{
 public:
   
   TT06_RRG_Reaction(const int numPoints, int ttType, const TT06_RRG_ReactionParms &parms);
   std::string methodName() const {return "TT06_RRG";}

   void calc(double dt,
             const Managed<ArrayView<double>> Vm,
             const Managed<ArrayView<double>> iStim,
             Managed<ArrayView<double>> dVm);
   void initializeMembraneVoltage(ArrayView<double> Vm);
   
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

   void initConst(const TT06_RRG_ReactionParms &parms); 
   double constants_[54]; 
   std::vector<TT06_RRG> cells_;
};

#endif
