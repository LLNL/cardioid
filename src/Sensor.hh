#ifndef SENSOR_HH
#define SENSOR_HH
#include <vector>

class Sensor
{
 public:
   virtual void print(double time, std::vector<double>& Vm) = 0;
   virtual void print(double time, std::vector<double>& Vm, std::vector<double>& dVm_r, std::vector<double>& dVm_d, std::vector<double>& dVm_e) = 0;
   virtual unsigned printRate(void) = 0;
   virtual bool printDerivs(void) = 0;
   virtual ~Sensor() {};
};

#endif
