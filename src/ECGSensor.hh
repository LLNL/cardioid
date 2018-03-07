#ifndef ECG_SENSOR_HH
#define ECG_SENSOR_HH

#include "Sensor.hh"
#include "VectorDouble32.hh"
#include "TransportCoordinator.hh"

struct ECGSensorParms
{
   int nFiles;
   int nSensorPoints;
   int stencilSize;
   std::string filename;
};

class Simulate;

class ECGSensor : public Sensor
{
 public:
   ECGSensor(const SensorParms& sp,
             const ECGSensorParms& p,
             const Simulate& sim);
   ~ECGSensor() {};

 private:

   void print(double time, int loop);
   void eval(double time, int loop);

   int nFiles_;
   unsigned nSensorPoints_;
   int stencilSize_;
   int nEval_;
   int dataOffset_;

   const TransportCoordinator<PinnedVector<double>>& VmTransport_;
   
   std::string filename_;
   std::vector<float> weight_;
   std::vector<int> VmOffset_;
   std::vector<float> data_;
};

#endif
