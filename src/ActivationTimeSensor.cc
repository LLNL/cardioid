#include "ActivationTimeSensor.hh"

#include "Anatomy.hh"
#include "pio.h"
#include "ioUtils.h"
#include "Simulate.hh"

#include <sstream>
#include <iomanip>

using namespace std;

ActivationTimeSensor::ActivationTimeSensor(const SensorParms& sp,
                                           const ActivationTimeSensorParms& p,
                                           const Anatomy& anatomy,
                                           const PotentialData& vdata)
: Sensor(sp),
  nLocal_(anatomy.nLocal()),
  nx_(anatomy.nx()),
  ny_(anatomy.ny()),
  nz_(anatomy.nz()),
  dx_(anatomy.dx()),
  dy_(anatomy.dy()),
  dz_(anatomy.dz()),
  vdata_(vdata),
  filename_(p.filename)
{
   activationTime_.resize(nLocal_, 0.0);
   activated_.resize(nLocal_, false);
   clear();
   cells_.reserve(nLocal_);
   for (unsigned ii=0; ii<nLocal_; ++ii)
      cells_.push_back(anatomy.globalTuple(ii));
}


void ActivationTimeSensor::print(double time, int loop)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   stringstream name;
   name << "snapshot."<<setfill('0')<<setw(12)<<loop;
   string fullname = name.str();
   if (myRank == 0)
      DirTestCreate(fullname.c_str());
   fullname += "/" + filename_;

   Long64 nGlobal;
   Long64 nLocal=nLocal_;
   MPI_Allreduce(&nLocal, &nGlobal, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
   
   PFILE* file = Popen(fullname.c_str(), "w", MPI_COMM_WORLD);

   char fmt[] = "%5d %5d %5d %18.12f";
   int lrec = 40;
   int nfields = 4; 


   if (myRank == 0)
   {
      // write header
      int nfiles;
      Pget(file,"ngroup",&nfiles);
      Pprintf(file, "activationTime FILEHEADER {\n");
      Pprintf(file, "  lrec = %d;\n", lrec);
      Pprintf(file, "  datatype = FIXRECORDASCII;\n");
      Pprintf(file, "  nrecords = %llu;\n", nGlobal);
      Pprintf(file, "  nfields = %d;\n", nfields);
      Pprintf(file, "  field_names = rx ry rz tActiv;\n");
      Pprintf(file, "  field_types = u u u f;\n" );
      Pprintf(file, "  nfiles = %u;\n", nfiles);
      Pprintf(file, "  dx = %f; dy = %f; dz = %f;\n", dx_, dy_, dz_);
      Pprintf(file, "  printRate = %d; evalRate = %d;\n", printRate(), evalRate());
      Pprintf(file, "  h = %4u  0    0\n", nx_);
      Pprintf(file, "        0    %4u  0\n", ny_);
      Pprintf(file, "        0    0    %4u;\n", nz_);
      Pprintf(file, "}\n\n");
   }

   int halfNx = nx_/2;
   int halfNy = ny_/2;
   int halfNz = nz_/2;
   char line[lrec+1];
   for (unsigned ii=0; ii<nLocal_; ++ii)
   {
      int ix = cells_[ii].x() - halfNx;
      int iy = cells_[ii].y() - halfNy;
      int iz = cells_[ii].z() - halfNz;

      int l = snprintf(line, lrec, fmt,
                       ix, iy, iz, activationTime_[ii]);
      
      for (; l < lrec - 1; l++) line[l] = (char)' ';
      line[l++] = (char)'\n';
      assert (l==lrec);
      Pwrite(line, lrec, 1, file);
   }
   Pclose(file);
}

void ActivationTimeSensor::eval(double time, int loop)
{
   for (unsigned ii=0; ii<nLocal_; ++ii)
   {
      if (activated_[ii]) continue;
      if ( (*vdata_.VmArray_)[ii] > 0 )
      {
         activated_[ii] = true;
         activationTime_[ii] = time;
      }
   }
}

void ActivationTimeSensor::clear()
{
   for (unsigned ii=0; ii<nLocal_; ++ii)
   {
      activated_[ii] = false;
      activationTime_[ii] = 0.0;
   }
}
