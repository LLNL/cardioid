#ifndef ACTIVATION_AND_RECOVERY_SENSOR
#define ACTIVATION_AND_RECOVERY_SENSOR

#include "Sensor.hh"
#include "Long64.hh"
#include <string>
#include <vector>

#include "Vector.hh"

class Anatomy;
class PotentialData;

struct ActivationAndRecoverySensorParms
{
   unsigned nFiles;
   std::string filename;
   double threshhold;
};

class ActivationAndRecoverySensor : public Sensor
{
 public:
   ActivationAndRecoverySensor(const SensorParms& sp,
                               const ActivationAndRecoverySensorParms& p,
                               const Anatomy& anatomy,
                               const PotentialData& vdata);
   ~ActivationAndRecoverySensor(){}
   
   void print(double time, int loop);
   void eval(double time, int loop);
   
 private:

   void clear();
   
   std::string filename_;

   unsigned nLocal_;
   unsigned nFiles_;
   double   threshhold_;

   std::vector<bool>                 active_;
   std::vector<std::vector<double> > activationTime_;
   std::vector<std::vector<double> > recoveryTime_;
   std::vector<Vector>               coords_;

   const PotentialData& vdata_;
};

#endif
