#ifndef DATAVORONOICOARSENING_HH
#define DATAVORONOICOARSENING_HH

#include "VoronoiCoarsening.hh"
#include "Sensor.hh"
#include "VectorDouble32.hh"

class PotentialData;
class Anatomy;

class DataVoronoiCoarsening : public Sensor
{
 private:
   VoronoiCoarsening coarsening_;
 
   std::string filename_;
   unsigned nFiles_;
   
   const Anatomy& anatomy_;
   const PotentialData& vdata_;

   MPI_Comm comm_;

   LocalSums avg_valcolors_;

   // eval times
   std::vector<double> times_;
   
   // average for each local color
   std::map<int,std::vector<float> > averages_;
   
   void computeColorAverages(const VectorDouble32& val);
   void writeAverages(const std::string& filename,
                      const double current_time,
                      const int current_loop)const;

 public:
   DataVoronoiCoarsening(const SensorParms& sp,
                         std::string filename,
                         unsigned nFiles,
                         const Anatomy& anatomy,
                         const std::vector<Long64>& gid,
                         const PotentialData& vdata,
                         const CommTable* commtable,
                         const double max_distance);
   void eval(double time, int loop);
   void print(double time, int loop);
};

#endif
