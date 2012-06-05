#ifndef ACTIVATION_TIME_SENSOR
#define ACTIVATION_TIME_SENSOR

#include "Sensor.hh"
#include <string>
#include <vector>

#include "Tuple.hh"

class Anatomy;
class PotentialData;

struct ActivationTimeSensorParms
{
   std::string filename;
};

class ActivationTimeSensor : public Sensor
{
 public:
   ActivationTimeSensor(const SensorParms& sp,
                        const ActivationTimeSensorParms& p,
                        const Anatomy& anatomy,
                        const PotentialData& vdata);
   ~ActivationTimeSensor(){}

   void print(double time, int loop);
   void eval(double time, int loop);
   

 private:

   void clear();
   
   std::string filename_;

   unsigned nLocal_;
   int nx_;
   int ny_;
   int nz_;
   double dx_;
   double dy_;
   double dz_;
   

   std::vector<bool>   activated_;
   std::vector<double> activationTime_;
   std::vector<Tuple>  cells_;

   const PotentialData& vdata_;
};

#endif
