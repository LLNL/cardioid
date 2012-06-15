#include "MaxDVSensor.hh"
#include "Simulate.hh"

#include <iomanip>
#include <ostream>
#include <fstream>

MaxDVSensor::MaxDVSensor(const SensorParms& sp,
                         const Anatomy& anatomy,
                         const PotentialData& vdata,
                         MPI_Comm comm, 
                         ostream* os)
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
                         string& filename)
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
      if( myRank_==0 )(*os_) << "#   Loop     Time         dVm"<<endl;
      first_time=false;
   }
   
   double maxdVdt=0.;
   double absMaxdVdt=0.;
   for (unsigned ii=0; ii<nlocal_; ++ii)
   {
      double dVdt=(*vdata_.dVmReaction_)[ii]+(*vdata_.dVmDiffusion_)[ii];
      double absdVdt=fabs( dVdt );
      if( absdVdt>absMaxdVdt )
      {
         absMaxdVdt=absdVdt;
         maxdVdt=dVdt;
      }
   }
   
   double maxMaxdVdt=0.;
   MPI_Reduce(&maxdVdt, &maxMaxdVdt, 1, MPI_DOUBLE, MPI_SUM, 0, comm_);

   if( myRank_==0 ){
      (*os_) << setw (8) << loop ;
      (*os_) << setw (9) << setprecision(3);
      os_->setf(ios_base::fixed,ios_base::floatfield);
      (*os_) << time ;
      (*os_) << setw (22)<< setprecision(15);
      (*os_) << maxMaxdVdt << endl;
   }
}
