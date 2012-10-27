#include "ECGSensor.hh"

#include <sstream>
#include <iomanip>
#include <cmath>

#include "pio.h"
#include "ioUtils.h"
#include "Simulate.hh"
#include "PerformanceTimers.hh"

using namespace std;
using PerformanceTimers::sensorEvalTimer;

ECGSensor::ECGSensor(const SensorParms& sp,
                     const ECGSensorParms& p,
                     const Simulate& sim)
: Sensor(sp),
  nFiles_(p.nFiles),
  nSensorPoints_(p.nSensorPoints),
  stencilSize_(p.stencilSize),
  filename_(p.filename),
  nEval_(0),
  Vm_(sim.vdata_.VmArray_)
{
   unsigned nLocal = sim.anatomy_.nLocal();
   nSensorPoints_ = min(nSensorPoints_, nLocal);
   
   weight_.resize(nSensorPoints_ * stencilSize_ * 3, 1.0);
   VmOffset_.resize(weight_.size());
   for (unsigned ii=0; ii<VmOffset_.size(); ++ii)
      VmOffset_[ii] = drand48()*nLocal;
   
   unsigned spaceNeeded = printRate()/evalRate() + 2;
   spaceNeeded *= 3;
   dataOffset_ = spaceNeeded;
   spaceNeeded *= nSensorPoints_;
   data_.resize(spaceNeeded);
}

void ECGSensor::print(double time, int loop)
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
 
   int nRec = nEval_+2; // Two dummy records for each sensor point.
   int lRec = 3*sizeof(float);
   for (unsigned ii=0; ii<nSensorPoints_; ++ii)
   {
      float* dataStartPtr = &data_[0] + ii*dataOffset_;
      Pwrite(dataStartPtr, lRec, nRec, file);
   }

   Pclose(file);
   nEval_ = 0;
}

void ECGSensor::eval(double time, int loop)
{
   startTimer(sensorEvalTimer);
   int index = 0;
   int dataIndex = (nEval_+2)*3;
   for (unsigned ii=0; ii<nSensorPoints_; ++ii)
   {
      for (unsigned jj=0; jj<3; ++jj)
      {
         float sum = 0.0;
         for (unsigned kk=0; kk<stencilSize_; ++kk)
         {
            sum += weight_[index] * Vm_[VmOffset_[index]];
            ++index;
         }
         data_[dataIndex+jj] = sum;
      }
      dataIndex += dataOffset_;
   }
   ++nEval_;
   stopTimer(sensorEvalTimer);
}

