#include "ECGSensor.hh"

#include <sstream>
#include <iomanip>
#include <cmath>

#include "pio.h"
#include "ioUtils.h"
#include "Simulate.hh"
#include "PerformanceTimers.hh"
#include "PioHeaderData.hh"
#include <cuda.h>
#include <cuda_runtime_api.h>

#define CUDA_VERIFY(x) do { cudaError_t error = x; if (error != cudaSuccess) { cout << error << endl; assert(error == cudaSuccess && #x ); } } while(0)

using namespace std;
using PerformanceTimers::sensorEvalTimer;

ECGSensor::ECGSensor(const SensorParms& sp,
                     const ECGSensorParms& p,
                     const Simulate& sim)
: Sensor(sp),
  nFiles_(p.nFiles),
  nSensorPoints_(p.nSensorPoints),
  stencilSize_(p.stencilSize),
  ecgNames(p.ecgNames),
  filename_(p.filename),
  nEval_(0),
  dVmDiffusionTransport_(sim.vdata_.dVmDiffusionTransport_)
{
    kECG=p.kconst;
    const int dim=3;
    nEcgPoints=p.ecgPoints.size()/dim;
    std::vector<double> ecgPoints(p.ecgPoints);
    ecgPointTransport_.resize(ecgPoints.size());
    auto ecgPointAccess = ecgPointTransport_.writeonly(CPU);
    copy(ecgPoints.begin(), ecgPoints.end(), ecgPointAccess.begin());
    calcInvR(sim);
   
    std::vector<double> ecgs(nEcgPoints, 0);
    ecgsTransport_.resize(ecgs.size());
    auto ecgAccess = ecgsTransport_.writeonly(CPU);
    copy(ecgs.begin(), ecgs.end(), ecgAccess.begin());
}

void ECGSensor::calcInvR(const Simulate& sim)
{
    Anatomy anatomy=sim.anatomy_;
    unsigned nlocal=anatomy.nLocal();
    int nx=anatomy.nx();
    int ny=anatomy.ny();
    int nz=anatomy.nz();

    double dx=anatomy.dx();
    double dy=anatomy.dy();
    double dz=anatomy.dz();

    kECG=kECG*dx*dy*dz;

    lazy_array<Long64>  gidsTransport_;
    gidsTransport_.resize(nlocal);
    auto gridAccess = gidsTransport_.writeonly(CPU);
    for(unsigned ii=0; ii<nlocal; ++ii){
        gridAccess[ii]=anatomy.gid(ii);
    }
    
    std::vector<double> invr(nlocal*nSensorPoints_,0.0);
    invrTransport_.resize(nlocal*nSensorPoints_);
    auto invrAccess=invrTransport_.writeonly(GPU);
    CUDA_VERIFY(cudaMemset(invrAccess.raw(), 0, sizeof(double)*invrAccess.size()));
    
    calcInvrCUDA(invrTransport_,
                 gidsTransport_,
                 ecgPointTransport_,
                 nEcgPoints,
                 nx, ny, nz,
                 dx, dy, dz);
        
 if(0){ // DEBUG to print out r with PIO
  int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   std::string filename_r="rData";
   int loop=0;
   double time=0;
   stringstream name;
   name << "snapshot."<<setfill('0')<<setw(12)<<loop;
   string fullname = name.str();
   if (myRank == 0)
      DirTestCreate(fullname.c_str());
   fullname += "/" + filename_r;

   int nFiles_r=4;
   PFILE* file = Popen(fullname.c_str(), "w", MPI_COMM_WORLD);
   if (nFiles_r > 0)
     PioSet(file, "ngroup", nFiles_r);

   char fmt[] = "%18.8g";
   int nFields = nEcgPoints;
   int lRec = 20*nFields;

   std::string fieldNames = "";
   std::string fieldTypes = "";
   std::string fieldUnits = "";

   for(unsigned ii=0; ii<nEcgPoints; ++ii)
   {
      fieldNames = fieldNames + "  " + ecgNames[ii];
      fieldTypes = fieldTypes + "  f";
      fieldUnits = fieldUnits + "  mv";
   }


   PioHeaderData header;
   header.objectName_ = "rData";
   header.className_  = "FILEHEADER";
   header.dataType_   = PioHeaderData::ASCII;
   header.nRecords_   = anatomy.nGlobal();
   header.lRec_       = lRec;
   header.nFields_    = nFields;
   header.fieldNames_ = fieldNames;
   header.fieldTypes_ = fieldTypes;
   header.fieldUnits_ = fieldUnits;

   if (myRank == 0) {
      header.writeHeader(file, loop, time);

   }
   auto invrT=invrTransport_.readonly(CPU);

      char line[lRec+1];

      for (unsigned ii=0; ii<nlocal; ++ii)
      {
         int l = 0;

         for(unsigned jj=0; jj<nEcgPoints; ++jj){
            int index = ii*nEcgPoints+jj;
            int ll = snprintf(line+l, lRec, fmt, 1.0/invrT[index]);
            l=l+ll;
         }

         for (; l < lRec - 1; l++) line[l] = (char)' ';
         line[l++] = (char)'\n';
         assert (l==lRec);
         Pwrite(line, lRec, 1, file);
      }


   Pclose(file);
 }
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

   char fmt[] = "%20.8g";
   int nFields = nEcgPoints+1;
   int lRec = 20*nFields;

   std::string fieldNames = "Loop";
   std::string fieldTypes = "d";
   std::string fieldUnits = "1";

   for(unsigned ii=0; ii<nEcgPoints; ++ii)
   {
      fieldNames = fieldNames + "  " + ecgNames[ii];
      fieldTypes = fieldTypes + "  f";
      fieldUnits = fieldUnits + "  mv";
   }
    

   PioHeaderData header;
   header.objectName_ = "ecgData";
   header.className_  = "FILEHEADER";
   header.dataType_   = PioHeaderData::ASCII;
   header.nRecords_   = saveLoops.size();
   header.lRec_       = lRec;
   header.nFields_    = nFields;
   header.fieldNames_ = fieldNames;
   header.fieldTypes_ = fieldTypes;
   header.fieldUnits_ = fieldUnits;
   header.addItem("printRate", printRate());
   header.addItem("evalRate", evalRate());

   if (myRank == 0) {
      header.writeHeader(file, loop, time);

      char line[lRec+1];
      
      for (unsigned ii=0; ii<saveLoops.size(); ++ii)
      {
         int l = snprintf(line, lRec, "%10d ", saveLoops[ii]);

         for(unsigned jj=0; jj<nEcgPoints; ++jj){
            int index = ii*nEcgPoints+jj;
            assert(index<saveEcgs.size());
	    int ll = snprintf(line+l, lRec, fmt, saveEcgs[index]);
            l=l+ll;
	 }
         //printf("l=%d, lRec=%d, ii=%d, ecgs=%f\n", l, lRec, ii, ecgs[ii]);
         for (; l < lRec - 1; l++) line[l] = (char)' ';
         line[l++] = (char)'\n';
         assert (l==lRec);
         Pwrite(line, lRec, 1, file);
      }

      // Clear up the loop and ecgs save values in vector
      saveLoops.clear();
      saveEcgs.clear();
      
   }

   Pclose(file);
}

void calcEcg(rw_mgarray_ptr<double> _ecgs,
             ro_mgarray_ptr<double> _invr,
             ro_mgarray_ptr<double> _dVmDiffusion,
             const int nEcgPoints) {
   auto ecgs = _ecgs.useOn(CPU);
   auto invr = _invr.useOn(CPU);
   auto dVmDiffusion = _dVmDiffusion.useOn(CPU);
   int nData=dVmDiffusion.size();
   for(unsigned ii=0; ii< nData; ii++){
      for (unsigned jj=0; jj< nEcgPoints; jj++){
         unsigned index=ii*nEcgPoints+jj;
         ecgs[jj]+=invr[index]*dVmDiffusion[ii];
      }
   }
}

void ECGSensor::eval(double time, int loop)
{
   startTimer(sensorEvalTimer);
   {   // zero out
   	auto ecgs = ecgsTransport_.writeonly(GPU);
        CUDA_VERIFY(cudaMemset(ecgs.raw(), 0, sizeof(double)*ecgs.size()));
   }
   
  if(0) { // DEBUG to print out invr and  dVmDiffusion
     if(loop<100 || loop%100==1){
        auto invr=invrTransport_.readonly(CPU);
        auto dVmDiffusion=dVmDiffusionTransport_.readonly(CPU);

   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   if(myRank==0){
       std::cout << "Loop "<< loop <<"\n\n\n";
	for(int ii=0; ii<100; ii++){
	    std::cout << invr[ii] << " ";
	    if((ii+1)%10==0) std::cout << "\n";
	}
        std::cout << "\n\n\n";
        for(int ii=0; ii<100; ii++){
            std::cout << dVmDiffusion[ii] << " ";
            if((ii+1)%10==0) std::cout << "\n";
        }
        std::cout << "\n\n\n";

   dump_GPU_data(ecgsTransport_,
               invrTransport_,
               dVmDiffusionTransport_, 
               nEcgPoints);
    }      

   }
  }

 if(0){ // DEBUG use the CPU version
   calcEcg(ecgsTransport_,
               invrTransport_,
               dVmDiffusionTransport_,
               nEcgPoints);
 }
 else
 {
   calcEcgCUDA(ecgsTransport_,
               invrTransport_,
               dVmDiffusionTransport_, 
               nEcgPoints);
 }
   
   auto ecgs = ecgsTransport_.readonly(CPU);
   const double* ecgsSendBuf=ecgs.raw();
   double ecgsRecvBuf[nEcgPoints];
   // MPI_Allreduce(ecgsSendBuf, ecgsRecvBuf, nEcgPoints, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);  
   // Only the Rank 0 stores the total ecg values
   MPI_Reduce(ecgsSendBuf, ecgsRecvBuf, nEcgPoints, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);  
 
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   if(myRank==0){   
        saveLoops.push_back(loop);
   	for(int ii=0; ii<ecgs.size(); ii++)
   	{
                double ecgValue=ecgsRecvBuf[ii]*kECG;
		saveEcgs.push_back(ecgValue);
   	}
   }
 
   stopTimer(sensorEvalTimer);
}

