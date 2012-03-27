#ifndef POINT_LIST_SENSOR_HH
#define POINT_LIST_SENSOR_HH

#include "Sensor.hh"
#include <vector>
#include <string>
#include <fstream>
#include "Long64.hh"

using namespace std;

class Anatomy;

struct PointListSensorParms
{
  vector<Long64> pointList;
  double startTime;
  double endTime;
   string filename;
   string dirname;
  int printDerivs;
};

class PointListSensor : public Sensor
{
 public:
   PointListSensor(const SensorParms& sp, const PointListSensorParms& p, const Anatomy& anatomy);
   ~PointListSensor();

   void print(double time, int loop, const vector<double>& Vm, const vector<double>& dVm_r, const vector<double>& dVm_d);
   void eval(double time, int loop,
             const std::vector<double>& Vm, const std::vector<double>& dVm_r,
             const std::vector<double>& dVm_d)
   {} // no eval function.
    
 private:
   void print(double time, const vector<double>& Vm);
   void printDerivs(double time, const vector<double>& Vm, const vector<double>& dVm_r,
                    const vector<double>& dVm_d);

   vector<Long64> localCells_;  // grid gids owned by this task
   vector<unsigned> sensorind_;      // corresponding local array index 
   vector<ofstream*> fout_loc_;
   double startTime_;
   double endTime_;
   bool printDerivs_;
};

#endif
