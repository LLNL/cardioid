#ifndef POINT_LIST_SENSOR_HH
#define POINT_LIST_SENSOR_HH

#include "Sensor.hh"
#include <vector>
#include <string>
#include <fstream>
#include "Long64.hh"

class Anatomy;
class PotentialData;

struct PointListSensorParms
{
   std::vector<Long64> pointList;
   std::string filename;
   std::string dirname;
   int printDerivs;
};

class PointListSensor : public Sensor
{
 public:
   PointListSensor(const SensorParms& sp, const PointListSensorParms& p, const Anatomy& anatomy, const PotentialData& vdata);
   ~PointListSensor();

   void print(double time, int loop);
   void eval(double time, int loop)
   {} // no eval function.
    
 private:
   void print(double time);
   void printDerivs(double time);

   std::vector<Long64> localCells_;  // grid gids owned by this task
   std::vector<unsigned> sensorind_;      // corresponding local array index 
   std::vector<std::ofstream*> fout_loc_;
   bool printDerivs_;
   
   const PotentialData& vdata_;
};

#endif
