#ifndef CA_SENSOR_HH
#define CA_SENSOR_HH

#include "VoronoiCoarsening.hh"
#include "Sensor.hh"

#include <string>
using namespace std;

class Reaction;
class Anatomy;

class CaAverageSensor : public Sensor
{
private:
   VoronoiCoarsening coarsening_;

   std::string filename_;

   const Reaction& reaction_;

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
   
   void computeColorAverages(const VectorDouble32& val);
   void writeAverages(const string& filename,
                      const double current_time,
                      const int current_loop)const;
public:
   CaAverageSensor(const SensorParms& sp,
            string filename,
            const Anatomy& anatomy,
            const vector<Long64>& gid,
            const Reaction& reaction,
            const CommTable* commtable);
   ~CaAverageSensor(){};
   
   void print(double time, int loop);
   void eval(double time, int loop);
   void bufferReactionData(const int loop);
   void bufferReactionData(const int begin, 
                           const int end, const int loop);
};


#endif
