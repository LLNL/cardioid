#ifndef SENSOR_HH
#define SENSOR_HH
#include <vector>

class Sensor
{
 public:
   virtual void print(double time, std::vector<double>& Vm) = 0;
   virtual unsigned printRate(void) = 0;
   virtual ~Sensor() {};
};

#endif
