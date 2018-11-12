#ifndef ECG_SENSOR_HH
#define ECG_SENSOR_HH

#include "Long64.hh"
#include "Sensor.hh"
#include "VectorDouble32.hh"
#include "lazy_array.hh"


void calcInvrCUDA(wo_mgarray_ptr<double> invr,
                  ro_mgarray_ptr<Long64> gids,
                  ro_mgarray_ptr<double> ecgPoints,
                  const int nEcgPoints,
                  const int nx, const int ny, const int nz,
                  const double dx, const double dy, const double dz);

void calcEcgCUDA(wo_mgarray_ptr<double> ecgs,
                 ro_mgarray_ptr<double> invr,
                 ro_mgarray_ptr<double> Vm,
                 const int nEcgPoints);
void dump_GPU_data(wo_mgarray_ptr<double> ecgs,
                   ro_mgarray_ptr<double> invr,
                   ro_mgarray_ptr<double> Vm,
                   const int nEcgPoints);

struct ECGSensorParms
{
   int nFiles;
   int nSensorPoints;
   int stencilSize;
   double kconst;
   std::string filename;
   std::vector<double> ecgPoints;
   std::vector<std::string> ecgNames;
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
   
   void calcInvR(const Simulate& sim);

   std::string filename_;
   
   int nFiles_;
   unsigned nSensorPoints_;
   int stencilSize_;
   int nEval_;
   int dataOffset_;
   std::vector<std::string> ecgNames;
   
   int nEcgPoints;

   double kECG;

   std::vector<int> saveLoops;
   std::vector<double> saveEcgs;

   const lazy_array<double>& dVmDiffusionTransport_;
   
   lazy_array<double> ecgPointTransport_;
   lazy_array<double> ecgsTransport_;
   lazy_array<double> invrTransport_;

};

#endif
