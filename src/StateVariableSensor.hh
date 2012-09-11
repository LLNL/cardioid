#ifndef STATE_VARIABLE_SENSOR_HH
#define STATE_VARIABLE_SENSOR_HH

#include "Sensor.hh"
#include <vector>
#include <string>
#include <fstream>
#include "Long64.hh"

using namespace std;

class Simulate;
class PotentialData;

struct StateVariableSensorParms
{
    Long64 gidCenter;
    double radius;
    vector<string> fieldList;
    double startTime;
    double endTime;
    string dirname;
    int printDerivs;
};

class StateVariableSensor : public Sensor
{
 public:
   StateVariableSensor(const SensorParms& sp, const StateVariableSensorParms& p, const Simulate& sim);
   ~StateVariableSensor();

   void print(double time, int loop);
   void eval(double time, int loop)
   {} // no eval function.
    
 private:
    void print(double time);

    const Simulate& sim_;
    Long64 gidCenter_;
    double radius_;
    vector<string> fieldNames_;
    vector<int> handles_;
    vector<Long64> localCells_;  // grid gids owned by this task
   vector<unsigned> sensorind_;      // corresponding local array index 
    vector<ofstream*> fout_loc_;
    double startTime_;
    double endTime_;
};

#endif
