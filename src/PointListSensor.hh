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
  unsigned printRate;
  int printDerivs;
};

class PointListSensor : public Sensor
{
 public:
   PointListSensor(const PointListSensorParms& p, const Anatomy& anatomy);
   ~PointListSensor();
   void print(double time, vector<double>& Vm);
   void print(double time, vector<double>& Vm, vector<double>& dVm_r, vector<double>& dVm_d, vector<double>& dVm_e);
   unsigned printRate(void) { return printRate_; };
   bool printDerivs(void) { return printDerivs_; };
    
 private:
   vector<unsigned> pointlist_loc_;  // grid gids owned by this task
   vector<unsigned> sensorind_;      // corresponding local array index 
   vector<ofstream*> fout_loc_;
   double startTime_;
   double endTime_;
   string filebase_;
   unsigned printRate_;
    bool printDerivs_;
};

#endif
