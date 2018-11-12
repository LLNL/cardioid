#include "MaxDVSensor.hh"
#include "Simulate.hh"

#include <iomanip>
#include <ostream>
#include <fstream>

using namespace std;

MaxDVSensor::MaxDVSensor(const SensorParms& sp,
                         const Anatomy& anatomy,
                         const PotentialData& vdata,
                         MPI_Comm comm, 
                         std::ostream* os)
:  Sensor(sp),
   vdata_(vdata),
   comm_(comm),
   os_(os)
{
   MPI_Comm comm_ = MPI_COMM_WORLD;
   MPI_Comm_rank(comm_, &myRank_);

   nlocal_=anatomy.nLocal();
   
   opened_file_=false;

   assert( os_!=0 );
}


MaxDVSensor::MaxDVSensor(const SensorParms& sp,
                         const Anatomy& anatomy,
                         const PotentialData& vdata,
                         MPI_Comm comm, 
                         std::string& filename)
:  Sensor(sp),
   vdata_(vdata),
   comm_(comm)
{
   MPI_Comm comm_ = MPI_COMM_WORLD;
   MPI_Comm_rank(comm_, &myRank_);

   nlocal_=anatomy.nLocal();
   
   if( filename == "cout" ){
      os_= &std::cout;
   }else{
      os_ = new ofstream(filename.data(),ios_base::out);
      opened_file_=true;
   }

   assert( os_!=0 );
}

void MaxDVSensor::print(double time, int loop)
{
   static bool first_time=true;
   
   if( first_time ){
      if( myRank_==0 )(*os_) << "#   Loop     Time         min. dVm        max. dVm"<<endl;
      first_time=false;
   }
   
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

   if( myRank_==0 ){
      (*os_) << setw (8) << loop ;
      (*os_) << setw (9) << setprecision(3);
      os_->setf(ios_base::fixed,ios_base::floatfield);
      (*os_) << time ;
      (*os_) << setw (22)<< setprecision(15)
             << minMindVdt;
      (*os_) << setw (22)<< setprecision(15)
             << maxMaxdVdt << endl;
   }
}
