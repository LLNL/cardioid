#ifndef OHaraRudy_HH
#define OHaraRudy_HH

#include <string>
#include <vector>
#include <map>
#include "CheckpointVarInfo.hh"
#include "OHaraRudy.h" 


static std::string OHaraRudyCurrentNames[]={"INaCai","INaCass","INaK","INab","ICab","IKb","IpCa","INaFast","INaL","Ito","IKr","IKs","IK1","ICa",""};
struct OHaraRudy_Parms
{
      std::vector<std::string> currentNames;
      std::vector<std::string> currentModels;
};
class OHaraRudy
{
 public:

   // There is a map of strings to VarHandle in getHandleMap that must
   // be maintained in parallel with this enum.  The value undefinedName
   // must be -1 and nVars must be last in this enum.

   OHaraRudy(int cellType, OHaraRudy_Parms &parms);
   virtual double calc(double dt, double Vm, double iStim);
   double defaultVoltage();
   static void getCheckpointInfo(std::vector<std::string>& fieldNames, std::vector<std::string>& fieldUnits);

   static int getVarHandle(const std::string& varName);
   void setValue(int varHandle, double value);

   double getValue(int handle) const;
   void getValue(const std::vector<int>& handle, std::vector<double>& value) const;
   static const std::string getUnit(const std::string& varName);
   virtual ~OHaraRudy(){}

 protected:

   static COMPONENTINFO *info_; 
   static HandleMap  handleMap_; 
   static int *privateStateOffset_; 
   static int *privateParmsOffset_; 
   static double *defaultStateENDO_;
   static double *defaultStateEPI_;
   static double *defaultStateM_;
   static double *defaultParmsENDO_;
   static double *defaultParmsEPI_;
   static double *defaultParmsM_;
   static HandleMap& getHandleMap();
   static int stateSize_; 
   static int parmsSize_; 
   void initConsts(int cellType,OHaraRudy_Parms &);
   void initState(int cellType);
   void initParms(int cellType);

   double defaultVoltage_;
   CELLPARMS *cellParms_;  
   double *state_;
};
class OHaraRudyDebug : public OHaraRudy
{
   public:
   OHaraRudyDebug(int cellType, OHaraRudy_Parms &parms) : OHaraRudy(cellType, parms){};
   virtual double calc(double dt, double Vm, double iStim);
};

#endif

