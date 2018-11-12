#ifndef STATE_VARIABLE_SENSOR_HH
#define STATE_VARIABLE_SENSOR_HH

#include "Sensor.hh"
#include <vector>
#include <string>
#include <map>
#include "Long64.hh"
#include "PioHeaderData.hh"


class Simulate;
class PotentialData;
class Reaction;

struct StateVariableSensorParms
{
   bool binaryOutput;
   bool allCells;
   bool allFields;
   unsigned nFiles;
   double radius;
   std::string filename;
   std::string cellListFilename;
   std::vector<Long64> cells;
   std::vector<std::string> fieldList;
};

class StateVariableSensor : public Sensor
{
   typedef std::map<Long64, unsigned> MapType; // gid to local array index
 public:
   StateVariableSensor(const SensorParms& sp, const StateVariableSensorParms& p, const Simulate& sim);
   ~StateVariableSensor();
   
   void print(double time, int loop);
   void eval(double time, int loop) {}; // no eval function.
   
 private:

   double getSimValue(int iCell, int varHandle);
   
   bool binaryOutput_;
   const Simulate& sim_;
   std::string filename_;
   std::string headerProlog_;
   unsigned lRec_;
   MapType localCells_;
   std::vector<int> handles_;
   PioHeaderData header_;
   const char* gidFormat_;
   const char* varFormat_;
};

#endif
