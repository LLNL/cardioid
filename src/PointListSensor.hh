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
  int printRate;
};

class PointListSensor : public Sensor
{
 public:
   PointListSensor(const PointListSensorParms& p, const Anatomy& anatomy);
   ~PointListSensor();
   void print();
   
 private:
   std::vector<unsigned> pointlist_loc_;
   std::vector<std::ofstream*> fout_loc_;
   double startTime_;
   double endTime_;
   std::string filebase_;
   int printRate_;
};

#endif
