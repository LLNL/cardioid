#ifndef MAXDV_SENSOR_HH
#define MAXDV_SENSOR_HH

#include "Sensor.hh"

#include <vector>
#include <string>
#include <iostream>
#include <mpi.h>

class Anatomy;
class PotentialData;

class MaxDVSensor : public Sensor
{
 private:
   const PotentialData& vdata_;
   unsigned nlocal_;

   MPI_Comm comm_;

   void print(double time);
   int nLocal_;
   int myRank_;
   std::ostream* os_;
   bool opened_file_;

 public:
   MaxDVSensor(const SensorParms& sp, const Anatomy& anatomy, const PotentialData& vdata, 
               MPI_Comm comm, std::ostream* os= &std::cout);
   MaxDVSensor(const SensorParms& sp, const Anatomy& anatomy, const PotentialData& vdata, 
               MPI_Comm comm, std::string& filename);
   ~MaxDVSensor()
   {
      if( opened_file_ )delete os_;
   };

   void print(double time, int loop);
   void eval(double time, int loop)
   {} // no eval function.    
};

#endif
