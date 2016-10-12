#include "Reaction.hh"
#include "object.h"
#include <vector>
class Anatomy;

namespace SimpleOHaraRudyTT06 
{

struct State 
{
   //EDIT_STATE
   double nai;
   double nass;
   double ki;
   double kss;
   double cai;
   double cass;
   double cansr;
   double cajsr;

   double mTT2;
   double hTT2;
   double jTT2;

   double mL;
   double hL;
   double hLp;
   double a;
   double iF;
   double iS;
   double ap;
   double iFp;
   double iSp;
   double d;
   double ff;
   double fs;
   double fcaf;
   double fcas;
   double jca;
   double nca;
   double ffp;
   double fcafp;
   double xrf;
   double xrs;
   double xs1;
   double xs2;
   double xk1;
   double Jrelnp;
   double Jrelp;
   double CaMKt;
};

struct PerCellFlags 
{
   //EDIT_PERCELL_FLAGS
};

struct PerCellParameters
{
   //EDIT_PERCELL_PARAMETERS
};

class ThisReaction : public Reaction
{
 public:
   ThisReaction(const Anatomy& anatomy);
   std::string methodName() const;
   
   void calc(double dt,
             const VectorDouble32& Vm,
             const std::vector<double>& iStim,
             VectorDouble32& dVm);
   void initializeMembraneVoltage(VectorDouble32& Vm);
   virtual void getCheckpointInfo(std::vector<std::string>& fieldNames,
                                  std::vector<std::string>& fieldUnits) const;
   virtual int getVarHandle(const std::string& varName) const;
   virtual void setValue(int iCell, int varHandle, double value);
   virtual double getValue(int iCell, int varHandle) const;
   virtual const std::string getUnit(const std::string& varName) const;

 public:
   //constant flags
   //EDIT_FLAGS
   int celltype;

   //EDIT_PARAMETERS
   double GCaB;

   //per-cell flags
   std::vector<PerCellFlags> perCellFlags_;
   std::vector<PerCellParameters> perCellParameters_;

 private:
   unsigned nCells_;
   std::vector<State> state_;
};

}

namespace scanReaction 
{
   Reaction* scanSimpleOHaraRudyTT06(OBJECT* obj, const Anatomy& anatomy);
}
