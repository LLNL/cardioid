#ifndef OHaraRudy_HH
#define OHaraRudy_HH

#include <string>
#include <vector>
#include <map>
#include "CheckpointVarInfo.hh"
#include "OHaraRudy.h" 


class OHaraRudy
{
 public:

   // There is a map of strings to VarHandle in getHandleMap that must
   // be maintained in parallel with this enum.  The value undefinedName
   // must be -1 and nVars must be last in this enum.
#include "OHaraRudyEnum.hh"

   OHaraRudy(int cellType);
   double calc(double dt, double Vm, double iStim);
   double defaultVoltage();
   static void getCheckpointInfo(std::vector<std::string>& fieldNames,
                                 std::vector<std::string>& fieldUnits);

   static int getVarHandle(const std::string& varName);
   void setValue(int varHandle, double value);

   double getValue(int handle) const;
   void getValue(const std::vector<int>& handle,
                 std::vector<double>& value) const;
   static const std::string getUnit(const std::string& varName);


   
 private:

   static HandleMap& getHandleMap();
   void initConsts(int cellType);
   void initStates(int cellType);


   double defaultVoltage_;
   CELLPARMS *cellParms_;  
   STATE state_;
};

#endif

