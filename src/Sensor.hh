#ifndef SENSOR_HH
#define SENSOR_HH
#include <vector>

struct SensorParms
{
   int evalRate;
   int printRate;
};


class Sensor
{
 public:
   Sensor(const SensorParms& p)
   : evalRate_(p.evalRate),
     printRate_(p.printRate)
   {}
   virtual ~Sensor() {};

   virtual void print(double time, int loop,
                      const std::vector<double>& Vm, const std::vector<double>& dVm_r,
                      const std::vector<double>& dVm_d) = 0;
   virtual void eval(double time, int loop,
                     const std::vector<double>& Vm, const std::vector<double>& dVm_r,
                     const std::vector<double>& dVm_d) = 0;
   int printRate(){return printRate_;}
   int evalRate(){return evalRate_;}

 private:
   int evalRate_;
   int printRate_;
};

#endif
