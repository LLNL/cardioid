#ifndef DATAVORONOICOARSENING_HH
#define DATAVORONOICOARSENING_HH

#include "VoronoiCoarsening.hh"
#include "Sensor.hh"

class DataVoronoiCoarsening : public Sensor
{
 private:
   VoronoiCoarsening coarsening_;
 
   std::string filename_;
   const Anatomy& anatomy_;

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
                     MPI_Comm comm);
   void eval(double time, int loop,
             const std::vector<double>& Vm, const std::vector<double>&,
             const std::vector<double>&);
   void print(double time, int loop,
              const std::vector<double>& Vm, const std::vector<double>&,
              const std::vector<double>&);
};

#endif
