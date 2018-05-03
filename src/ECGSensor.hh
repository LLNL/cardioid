#ifndef ECG_SENSOR_HH
#define ECG_SENSOR_HH

#include "Long64.hh"
#include "Sensor.hh"
#include "VectorDouble32.hh"
#include "TransportCoordinator.hh"


void calcInvrCUDA(OnDevice<ArrayView<double>> invr,
                  OnDevice<ConstArrayView<Long64>> gids,
                  OnDevice<ConstArrayView<double>> ecgPoints,
                  const int nEcgPoints,
                  const int nx, const int ny, const int nz,
                  const double dx, const double dy, const double dz);

void calcEcgCUDA(OnDevice<ArrayView<double>> ecgs,
                 OnDevice<ConstArrayView<double>> invr,
                 OnDevice<ConstArrayView<double>> Vm,
                 const int nEcgPoints);

struct ECGSensorParms
{
   int nFiles;
   int nSensorPoints;
   int stencilSize;
   double kconst;
   std::string filename;
   PinnedVector<double> ecgPoints;
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
   
   int nEcgPoints;

   double kECG;

   std::vector<int> saveLoops;
   std::vector<double> saveEcgs;

   const TransportCoordinator<PinnedVector<double>>& VmTransport_;
   
   TransportCoordinator<PinnedVector<double> > ecgPointTransport_;

   TransportCoordinator<PinnedVector<double> > ecgsTransport_;
   
   TransportCoordinator<PinnedVector<double> > invrTransport_;

};

#endif
