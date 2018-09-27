#include "ActivationTimeSensor.hh"

#include "Anatomy.hh"
#include "PioHeaderData.hh"
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
  filename_(p.filename),
  nFiles_(p.nFiles)
{
   activationTime_.resize(nLocal_, 0.0);
   activated_.resize(nLocal_, false);
   clear();
   cells_.reserve(nLocal_);
   for (unsigned ii=0; ii<nLocal_; ++ii)
      cells_.push_back(anatomy.gid(ii));
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
   if (nFiles_ > 0)
     PioSet(file, "ngroup", nFiles_);
   
   char fmt[] = "%12llu %18.12f";
   int lRec = 32;
   int nFields = 2; 

   PioHeaderData header;
   header.objectName_ = "activationTime";
   header.className_  = "FILEHEADER";
   header.dataType_   = PioHeaderData::ASCII;
   header.nRecords_   = nGlobal;
   header.lRec_       = lRec;
   header.nFields_    = nFields;
   header.fieldNames_ = "gid tActiv";
   header.fieldTypes_ = "u f";
   header.fieldUnits_ = "1 ms";
   header.addItem("nx", nx_);
   header.addItem("ny", ny_);
   header.addItem("nz", nz_);
   header.addItem("dx", dx_);
   header.addItem("dy", dy_);
   header.addItem("dz", dz_);
   header.addItem("printRate", printRate());
   header.addItem("evalRate", evalRate());
   
   if (myRank == 0)
      header.writeHeader(file, loop, time);

   char line[lRec+1];
   for (unsigned ii=0; ii<nLocal_; ++ii)
   {
      int l = snprintf(line, lRec, fmt,
                       cells_[ii], activationTime_[ii]);
      
      for (; l < lRec - 1; l++) line[l] = (char)' ';
      line[l++] = (char)'\n';
      assert (l==lRec);
      Pwrite(line, lRec, 1, file);
   }
   Pclose(file);
}

void ActivationTimeSensor::eval(double time, int loop)
{
   ro_array_ptr<double> VmArray = vdata_.VmTransport_.useOn(CPU);
   for (unsigned ii=0; ii<nLocal_; ++ii)
   {
      if (activated_[ii]) continue;
      if (VmArray[ii] > 0 )
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
