#ifndef DATAVORONOICOARSENING_HH
#define DATAVORONOICOARSENING_HH

#include "VoronoiCoarsening.hh"
#include "Sensor.hh"

class PotentialData;
class Anatomy;

class DataVoronoiCoarsening : public Sensor
{
 private:
   VoronoiCoarsening coarsening_;
 
   std::string filename_;
   const Anatomy& anatomy_;
   const PotentialData& vdata_;

   MPI_Comm comm_;

   LocalSums avg_valcolors_;

   void computeColorAverages(const std::vector<double>& val);
   void writeAverages(const std::string& filename,
                      const double current_time,
                      const int current_loop)const;

 public:
   DataVoronoiCoarsening(const SensorParms& sp,
                     std::string filename,
                     const Anatomy& anatomy,
                     const std::vector<Long64>& gid,
                     const PotentialData& vdata,
                     MPI_Comm comm);
   void eval(double time, int loop);
   void print(double time, int loop);
};

#endif
