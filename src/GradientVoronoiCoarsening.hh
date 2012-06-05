#ifndef GRADIENTVORONOICOARSENING_HH
#define GRADIENTVORONOICOARSENING_HH

#include "VoronoiCoarsening.hh"
#include "Sensor.hh"

class Anatomy;
class PotentialData;

class GradientVoronoiCoarsening : public Sensor
{
 private:
   VoronoiCoarsening coarsening_;
 
   std::string filename_;
   const Anatomy& anatomy_;
   const PotentialData& vdata_;

   MPI_Comm comm_;

   std::vector<double> dx_;
   std::vector<double> dy_;
   std::vector<double> dz_;

   LocalSums valcolors_;
   LocalSums valMat00_;
   LocalSums valMat01_;
   LocalSums valMat02_;
   LocalSums valMat11_;
   LocalSums valMat12_;
   LocalSums valMat22_;
   LocalSums valRHS0_;
   LocalSums valRHS1_;
   LocalSums valRHS2_;

   void writeLeastSquareGradients(const std::string& filename,
                                  const double current_time,
                                  const int current_loop)const;
   void computeLSsystem(const std::vector<double>& val);
   void computeColorCenterValues(const std::vector<double>& val);

 public:
   GradientVoronoiCoarsening(const SensorParms& sp,
                     std::string filename,
                     const Anatomy& anatomy,
                     const std::vector<Long64>& gid,
                     const PotentialData& vdata_,
                     MPI_Comm comm);
   void eval(double time, int loop);
   void print(double time, int loop);
};

#endif
