#ifndef Grandi_REACTION_HH
#define Grandi_REACTION_HH

#include <string>
#include <vector>
#include <map>
#include "CheckpointVarInfo.hh"
#include "Grandi.h" 
#include "Reaction.hh"
class Anatomy;
class Grandi;
class BucketOfBits;
struct Grandi_Parms
{
   CELLTYPES cellType;
      std::vector<std::string> currentNames;
      std::vector<std::string> currentModels;
};

class Grandi_Reaction : public Reaction
{
 public:
   
   Grandi_Reaction(const int numPoints, Grandi_Parms &parms);
   std::string methodName() const {return "Grandi";}
   void calc(double dt,
             const VectorDouble32& Vm,
             const std::vector<double>& iStim,
             VectorDouble32& dVm);
   void initializeMembraneVoltage(VectorDouble32& Vm);

   /** Functions needed for checkpoint/restart */
   void getCheckpointInfo(std::vector<std::string>& fieldNames, std::vector<std::string>& fieldUnits) const;
   int getVarHandle(const std::string& varName) const;
   void setValue(int iCell, int varHandle, double value);
   double getValue(int iCell, int varHandle) const;
   void getValue(int iCell,
                 const std::vector<int>& handle,
                 std::vector<double>& value) const;
   const std::string getUnit(const std::string& varName) const;
   
 protected:
   
   static HandleMap  handleMap_;
 private:

   int indexS_;
   int nCells_; 
   CELLTYPES cellType;
   void makeHandleMap(); 
};

#endif
