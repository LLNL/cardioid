#ifndef DVTHRESH_SENSOR_HH
#define DVTHRESH_SENSOR_HH

#include "Sensor.hh"

#include <vector>
#include <string>
#include <iostream>
#include <mpi.h>

class Anatomy;
class PotentialData;

class DVThreshSensor : public Sensor
{
 private:
   const PotentialData& vdata_;
   MPI_Comm comm_;

   int nLocal_;
   int myRank_;
   int nlocal_;
   double threshold_;
    
 public:
   DVThreshSensor(const SensorParms& sp, const Anatomy& anatomy,const PotentialData& vdata,
                  MPI_Comm comm);
   ~DVThreshSensor() {};

   void print(double time, int loop) {}; // no print function
   void eval(double time, int loop);
};

#endif
