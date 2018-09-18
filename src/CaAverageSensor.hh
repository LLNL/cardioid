#ifndef CA_SENSOR_HH
#define CA_SENSOR_HH

#include "VoronoiCoarsening.hh"
#include "Sensor.hh"
#include "VectorDouble32.hh"

#include <string>

class ReactionManager;
class Anatomy;

class CaAverageSensor : public Sensor
{
private:
   VoronoiCoarsening coarsening_;

   std::string filename_;
   unsigned nFiles_;

   const ReactionManager& reaction_;

   VectorDouble32 buffer_val_;
   int loop_buffer_;
   
   unsigned nx_;
   unsigned ny_;
   unsigned nz_;
   unsigned nlocal_;
   
   int ca_handle_;

   MPI_Comm comm_;

   LocalSums avg_valcolors_;

   // eval times
   std::vector<double> times_;
   
   // average for each local color
   std::map<int,std::vector<float> > averages_;
   
   void computeColorAverages(ro_array_ptr<double> val);
   void writeAverages(const std::string& filename,
                      const double current_time,
                      const int current_loop)const;
public:
   CaAverageSensor(const SensorParms& sp,
                   std::string filename,
                   unsigned nFiles,
                   const Anatomy& anatomy,
                   std::vector<Long64>& gid,
                   const ReactionManager& reaction,
                   const CommTable* commtable);
   ~CaAverageSensor(){};
   
   void print(double time, int loop);
   void eval(double time, int loop);
   void bufferReactionData(const int loop);
   void bufferReactionData(const int begin, 
                           const int end, const int loop);
};


#endif
