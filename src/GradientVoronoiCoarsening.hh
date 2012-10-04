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
   
   double max_distance_;

   std::vector<double> dx_;
   std::vector<double> dy_;
   std::vector<double> dz_;
   
   std::vector<int> colored_cells_;

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
   void computeLSsystem(const VectorDouble32& val);
   void computeColorCenterValues(const VectorDouble32& val);

 public:
   GradientVoronoiCoarsening(const SensorParms& sp,
                     std::string filename,
                     const Anatomy& anatomy,
                     const std::vector<Long64>& gid,
                     const PotentialData& vdata_,
                     MPI_Comm comm,
                     const double max_distance);
   void eval(double time, int loop);
   void print(double time, int loop);
};

#endif
