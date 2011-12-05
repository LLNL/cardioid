#ifndef POINT_LIST_SENSOR_HH
#define POINT_LIST_SENSOR_HH

#include "Sensor.hh"
#include <vector>
#include <string>
#include <fstream>
using namespace std;

class Anatomy;

struct PointListSensorParms
{
  vector<unsigned> pointlist;
  double startTime;
  double endTime;
  string filebase;
  int printDerivs;
};

class PointListSensor : public Sensor
{
 public:
   PointListSensor(const SensorParms& sp, const PointListSensorParms& p, const Anatomy& anatomy);
   ~PointListSensor();

   void print(double time, int loop, const vector<double>& Vm, const vector<double>& dVm_r, const vector<double>& dVm_d, const vector<double>& dVm_e);
   void eval(double time, int loop,
             const std::vector<double>& Vm, const std::vector<double>& dVm_r,
             const std::vector<double>& dVm_d, const std::vector<double>& dVm_e)
   {} // no eval function.
    
 private:
   void print(double time, const vector<double>& Vm);
   void printDerivs(double time, const vector<double>& Vm, const vector<double>& dVm_r,
                    const vector<double>& dVm_d, const vector<double>& dVm_e);

   vector<unsigned> pointlist_loc_;  // grid gids owned by this task
   vector<unsigned> sensorind_;      // corresponding local array index 
   vector<ofstream*> fout_loc_;
   double startTime_;
   double endTime_;
   string filebase_;
   bool printDerivs_;
};

#endif
