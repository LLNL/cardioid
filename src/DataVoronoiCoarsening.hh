#ifndef DATAVORONOICOARSENING_HH
#define DATAVORONOICOARSENING_HH

#include "VoronoiCoarsening.hh"
#include "Sensor.hh"
#include "lazy_array.hh"

class PotentialData;
class Anatomy;

//Includes hacks from JPC to get Coarsened Activation Times out (changes marked with //AT-HACK), would eventually like to make it less hacky and interface with object.data (one bad thing I am doing is disabling the capability that was in previous version of DataVoronoiCoarsening code that allows you to print out data from multiple time points when evalrate<printrate, would be better if this was an object.data option.  I'm disabling this b/c to find Activation Times you'd want to evaluate often but print less frequently, which triggers printing out data from multiple time points in coarsened_anatomy#* files, which crashes the post-processing anatomy2ensight visualization tool as well as ECG code)

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
   
   void computeColorAverages(ro_array_ptr<double> val);
   void writeAverages(const std::string& filename,
                      const double current_time,
                      const int current_loop)const;
   void writeAveragesAT(const std::string& filename,
                      const double current_time,
                      const int current_loop)const;	//AT-HACK, this whole function is new, and is a sister function of writeAverages, but this one tells cardioid how to write out coarsened AT, instead of coarsened Vm.  For confusing parts I've added further AT-HACK comments within this function.	      

   std::vector<bool> active_;					//AT-HACK, active status of select gids in sensor.txt?, true or false
   std::vector<std::vector<double> > activationTime_;		//AT-HACK, activation time in ms (normalized to initation of simulation at t=0 ms) of select gids in sensor.txt
   void clear();						//AT-HACK, function that sets active to false and AT to -1000 ms for all select gids in sensor.txt

 public:
   DataVoronoiCoarsening(const SensorParms& sp,
                         std::string filename,
                         unsigned nFiles,
                         const Anatomy& anatomy,
                         std::vector<Long64>& gid,
                         const PotentialData& vdata,
                         const CommTable* commtable,
                         const double max_distance);
   void eval(double time, int loop);
   void print(double time, int loop);
};

#endif
