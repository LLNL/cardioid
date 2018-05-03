#include "ECGSensor.hh"

#include <sstream>
#include <iomanip>
#include <cmath>

#include "pio.h"
#include "ioUtils.h"
#include "Simulate.hh"
#include "PerformanceTimers.hh"
#include "PioHeaderData.hh"

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
  dVmDiffusionTransport_(sim.vdata_.dVmDiffusionTransport_)
{
    kECG=p.kconst;
    const int dim=3;
    nEcgPoints=p.ecgPoints.size()/dim;
    PinnedVector<double> ecgPoints(p.ecgPoints);
    ecgPointTransport_.setup(std::move(ecgPoints));
    calcInvR(sim);
   
    PinnedVector<double> ecgs(nEcgPoints, 0);
    ecgsTransport_.setup(std::move(ecgs));
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

    TransportCoordinator<PinnedVector<Long64> > gidsTransport_;
    PinnedVector<Long64> gids(nlocal, 0);
    for(unsigned ii=0; ii<nlocal; ++ii){
        gids[ii]=anatomy.gid(ii);
    }
    gidsTransport_.setup(std::move(gids));
    
    PinnedVector<double> invr(nlocal*nSensorPoints_,0.0);
    invrTransport_.setup(std::move(invr));
    
    calcInvrCUDA(invrTransport_,
                 gidsTransport_,
                 ecgPointTransport_,
                 nEcgPoints,
                 nx, ny, nz,
                 dx, dy, dz);
        
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

   char fmt[] = "%20.8f";
   int nFields = nEcgPoints+1;
   int lRec = 20*nFields;

   std::string fieldNames = "Loop";
   std::string fieldTypes = "d";
   std::string fieldUnits = "1";

   for(unsigned ii=0; ii<nEcgPoints; ++ii)
   {
	fieldNames = fieldNames +"   ecgPoint"+std::to_string(ii+1);
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

void ECGSensor::eval(double time, int loop)
{
   startTimer(sensorEvalTimer);
   {   // zero out
   	ArrayView<double> ecgs = ecgsTransport_;
   	for(int ii=0; ii<ecgs.size(); ii++)
   	{
        	ecgs[ii]=0.0;
   	}
   }
   
   calcEcgCUDA(ecgsTransport_,
               invrTransport_,
               dVmDiffusionTransport_, 
               nEcgPoints);
   
   ArrayView<double> ecgs = ecgsTransport_;
   double* ecgsSendBuf=&ecgs[0];
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

