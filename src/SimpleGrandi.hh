#include "Reaction.hh"
#include "object.h"
#include <vector>
class Anatomy;

namespace SimpleGrandi 
{

struct State 
{
   //EDIT_STATE
   double m;   
   double h;   
   double j;   
   double mL;   
   double hL;   
   double d;   
   double f;   
   double fcaBj;   
   double fcaBsl;   
   double xtof;   
   double ytof;   
   double xkr;   
   double xks;   
   double xkur;   
   double ykur;   
   double RyRr;   
   double RyRo;   
   double RyRi;   
   double NaBj;   
   double NaBsl;   
   double TnCL;   
   double TnCHc;   
   double TnCHm;   
   double CaM;   
   double Myc;   
   double Mym;   
   double SRB;   
   double SLLj;   
   double SLLsl;   
   double SLHj;   
   double SLHsl;   
   double Csqnb;   
   double Naj;   
   double Nasl;   
   double Nai;   
   double Ki;   
   double Casr;   
   double Caj;   
   double Casl;   
   double Cai;   
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
   bool AF;
   bool ISO;
   bool RA;

   //constant parameters
   //EDIT_PARAMETERS
   double ks;
   double Vmax_SRCaP;
   double GCaB;
   double GClCa;
   double GClB;
   double gkp;
   double IbarNaK;
   double IbarSLCaP;

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
   Reaction* scanSimpleGrandi(OBJECT* obj, const Anatomy& anatomy);
}

  
