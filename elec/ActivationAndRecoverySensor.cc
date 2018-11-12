#include "ActivationAndRecoverySensor.hh"

#include "Anatomy.hh"
#include "pio.h"
#include "ioUtils.h"
#include "Simulate.hh"

#include <sstream>
#include <iomanip>

using namespace std;

ActivationAndRecoverySensor::ActivationAndRecoverySensor(
   const SensorParms& sp,
   const ActivationAndRecoverySensorParms& p,
   const Anatomy& anatomy,
   const PotentialData& vdata)
: Sensor(sp),
  nLocal_(anatomy.nLocal()),
  vdata_(vdata),
  filename_(p.filename),
  nFiles_(p.nFiles),
  threshhold_(p.threshhold)
{
   activationTime_.resize(nLocal_);
   recoveryTime_.resize(nLocal_);
   active_.resize(nLocal_, false);
   clear();

   IndexToTuple indexToTuple(anatomy.nx(), anatomy.ny(), anatomy.nz());
   coords_.reserve(nLocal_);

   for (unsigned ii=0; ii<nLocal_; ++ii)
   {
      Tuple tmp = indexToTuple(anatomy.gid(ii));
      coords_.push_back(Vector(tmp.x()*anatomy.dx() + anatomy.dx()/2.0,
                               tmp.y()*anatomy.dy() + anatomy.dy()/2.0,
                               tmp.z()*anatomy.dz() + anatomy.dz()/2.0)); 
   }
}


void ActivationAndRecoverySensor::print(double time, int loop)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   
   stringstream name;
   name << "snapshot."<<setfill('0')<<setw(12)<<loop;
   string fullname = name.str();
   if (myRank == 0)
      DirTestCreate(fullname.c_str());
   fullname += "/" + filename_;

   PFILE* file = Popen(fullname.c_str(), "w", MPI_COMM_WORLD);
   if (nFiles_ > 0)
     PioSet(file, "ngroup", nFiles_);

   for (unsigned ii=0; ii<nLocal_; ++ii)
   {
      assert( activationTime_[ii].size() - recoveryTime_[ii].size() <= 1);
      Pprintf(file, "%6.4f %6.4f %6.4f", coords_[ii].x(), coords_[ii].y(), coords_[ii].z());
      for (unsigned jj=0; jj<activationTime_[ii].size(); ++jj)
      {
         Pprintf(file, " %8.3f", activationTime_[ii][jj]);
         if (recoveryTime_[ii].size() > jj)
            Pprintf(file, " %8.3f", recoveryTime_[ii][jj]);
      }
      Pprintf(file, "\n");
   }
   
   Pclose(file);
   clear();
}

void ActivationAndRecoverySensor::eval(double time, int loop)
{
   ro_array_ptr<double> VmArray = vdata_.VmTransport_.useOn(CPU);
   for (unsigned ii=0; ii<nLocal_; ++ii)
   {
      if (active_[ii] && VmArray[ii] < threshhold_ )
      {
         active_[ii] = false;
         recoveryTime_[ii].push_back(time);
      }
      if (! active_[ii] && VmArray[ii] > threshhold_ )
      {
         active_[ii] = true;
         activationTime_[ii].push_back(time);
         continue;
      }
   }
}

void ActivationAndRecoverySensor::clear()
{
   for (unsigned ii=0; ii<nLocal_; ++ii)
   {
      active_[ii] = false;
//       if (vdata_.VmArray_[ii] > threshhold_)
//          active_[ii] = true;
      activationTime_[ii].clear();
      activationTime_[ii].reserve(10);
      recoveryTime_[ii].clear();
      recoveryTime_[ii].reserve(10);
   }
}
