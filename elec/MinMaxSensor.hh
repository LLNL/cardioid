#ifndef MINMAX_SENSOR_HH
#define MINMAX_SENSOR_HH

#include "Sensor.hh"
#include <vector>
#include <string>
#include <fstream>
#include "Long64.hh"

using namespace std;

class Anatomy;
class PotentialData;

struct MinMaxSensorParms
{
   string filename;
   string dirname;
};

class MinMaxSensor : public Sensor
{
 public:
   MinMaxSensor(const SensorParms& sp, const MinMaxSensorParms& p, const Anatomy& anatomy, const PotentialData& vdata);
   ~MinMaxSensor();

   void print(double time, int loop);
   void eval(double time, int loop)
   {} // no eval function.
    
 private:
    void print(double time);
    int nLocal_;
    int myRank_;
    ofstream* fout_;
    bool printDerivs_;
   
    const PotentialData& vdata_;
};

#endif
