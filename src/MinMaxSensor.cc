#include "MinMaxSensor.hh"
#include "Anatomy.hh"
#include "ioUtils.h"
#include "Simulate.hh"

#include <mpi.h>
#include <cmath>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

MinMaxSensor::MinMaxSensor(const SensorParms& sp, const MinMaxSensorParms& p, 
                                 const Anatomy& anatomy,
                                 const PotentialData& vdata)
: Sensor(sp),
  vdata_(vdata),
  nLocal_(anatomy.nLocal())
{
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &myRank_);

  string filename = p.dirname + "/" + p.filename;
      
  if (myRank_ == 0)
  {
     DirTestCreate(p.dirname.c_str());
     fout_ = new ofstream;
     fout_->open(filename.c_str(),ofstream::out);
     fout_->setf(ios::scientific,ios::floatfield);
     (*fout_) << "#    time   min V_m    max V_m    max-min" << endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

MinMaxSensor::~MinMaxSensor()
{
   fout_->close();
   delete fout_;
}

void MinMaxSensor::print(double time, int /*loop*/)
{
   // find local min/max voltages
   ro_array_ptr<double> VmArray = vdata_.VmTransport_.useOn(CPU);
   double vmin_loc = VmArray[0];
   double vmax_loc = VmArray[0];
   for (unsigned ii=1; ii<nLocal_; ++ii)
   {
      if ( VmArray[ii] > vmax_loc )
         vmax_loc = VmArray[ii];
      if ( VmArray[ii] < vmin_loc )
         vmin_loc = VmArray[ii];
   }
   
   // MPI_Allreduce over all tasks to get global min/max
   double vmin, vmax;
   MPI_Allreduce(&vmin_loc, &vmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
   MPI_Allreduce(&vmax_loc, &vmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   
   if (myRank_ == 0)
   {
      (*fout_) << setprecision(10) << " " << time << "     " << vmin << "      " << vmax << "    " << vmax-vmin << endl;
      //ewd DEBUG cout << setprecision(10) << " " << time << "     " << vmin << "      " << vmax << endl;
   }
}
