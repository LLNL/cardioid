#ifndef ACTIVATION_TIME_SENSOR
#define ACTIVATION_TIME_SENSOR

#include "Sensor.hh"
#include <string>
#include <vector>

#include "Tuple.hh"

class Anatomy;

struct ActivationTimeSensorParms
{
   std::string filename;
};

class ActivationTimeSensor : public Sensor
{
 public:
   ActivationTimeSensor(const SensorParms& sp,
                        const ActivationTimeSensorParms& p,
                        const Anatomy& anatomy);
   ~ActivationTimeSensor(){}

   void print(double time, int loop,
              const std::vector<double>& Vm, const std::vector<double>& dVm_r,
              const std::vector<double>& dVm_d, const std::vector<double>& dVm_e);
   void eval(double time, int loop,
             const std::vector<double>& Vm, const std::vector<double>& dVm_r,
             const std::vector<double>& dVm_d, const std::vector<double>& dVm_e);
   

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
};

#endif
