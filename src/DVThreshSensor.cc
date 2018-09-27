#include "DVThreshSensor.hh"
#include "Simulate.hh"

DVThreshSensor::DVThreshSensor(const SensorParms& sp,
                         const Anatomy& anatomy,
                         const PotentialData& vdata,
                               MPI_Comm comm)
:  Sensor(sp),
   vdata_(vdata),
   comm_(comm)
{
   MPI_Comm comm_ = MPI_COMM_WORLD;
   MPI_Comm_rank(comm_, &myRank_);
   nlocal_=anatomy.nLocal();
   threshold_ = sp.value;
}

void DVThreshSensor::eval(double time, int loop)
{
   double maxdVdt=-10000.;
   double mindVdt= 10000.;
   ro_array_ptr<double> dVmDiffusion = vdata_.dVmDiffusionTransport_.useOn(CPU);
   ro_array_ptr<double> dVmReaction = vdata_.dVmReactionTransport_.useOn(CPU);

   for (unsigned ii=0; ii<nlocal_; ++ii)
   {
      double dVdt=dVmReaction[ii]+dVmDiffusion[ii];
      if( dVdt>maxdVdt )
      {
         maxdVdt=dVdt;
      }
      if( dVdt<mindVdt )
      {
         mindVdt=dVdt;
      }
   }
   
   double maxMaxdVdt=0.;
   MPI_Reduce(&maxdVdt, &maxMaxdVdt, 1, MPI_DOUBLE, MPI_MAX, 0, comm_);
   double minMindVdt=0.;
   MPI_Reduce(&mindVdt, &minMindVdt, 1, MPI_DOUBLE, MPI_MIN, 0, comm_);

   // if abs of max and min are both less than sp.value, stop simulation
   if (fabs(maxMaxdVdt) < threshold_ && fabs(minMindVdt) < threshold_)
   {
      if (myRank_ == 0)
         std::cout << "DVThreshold sensor:  maxdVdt = " << maxMaxdVdt << ", mindVdt = " << minMindVdt << ", threshold = " << threshold_ << " reached." << std::endl;
      exit(1);
   }
}
