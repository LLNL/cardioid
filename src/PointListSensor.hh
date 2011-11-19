#ifndef POINT_LIST_SENSOR_HH
#define POINT_LIST_SENSOR_HH

#include "Sensor.hh"
#include <vector>
#include <string>
#include <fstream>

class Anatomy;

struct PointListSensorParms
{
  std::vector<unsigned> pointlist;
  double startTime;
  double endTime;
  std::string filebase;
  unsigned printRate;
};

class PointListSensor : public Sensor
{
 public:
   PointListSensor(const PointListSensorParms& p, const Anatomy& anatomy);
   ~PointListSensor();
   void print(double time, std::vector<double>& Vm);
   unsigned printRate(void) { return printRate_; };
    
 private:
   std::vector<unsigned> pointlist_loc_;  // grid gids owned by this task
   std::vector<unsigned> sensorind_;      // corresponding local array index 
   std::vector<std::ofstream*> fout_loc_;
   double startTime_;
   double endTime_;
   std::string filebase_;
   unsigned printRate_;
};

#endif
