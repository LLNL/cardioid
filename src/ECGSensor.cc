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
  VmTransport_(sim.vdata_.VmTransport_)
{
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
 
   int nRec = nEval_+2; // Two dummy records for each sensor point.
   int lRec = 3*sizeof(float);

   ArrayView<double> ecgs = ecgsTransport_;
   for (unsigned ii=0; ii<ecgs.size(); ++ii)
   {
      double* dataStartPtr = &ecgs[0] + ii;
      Pwrite(dataStartPtr, lRec, nRec, file);
   }

   Pclose(file);
   nEval_ = 0;
}

void ECGSensor::eval(double time, int loop)
{
   ConstArrayView<double> Vm = VmTransport_;
   startTimer(sensorEvalTimer);
   
   calcEcgCUDA(ecgsTransport_,
               invrTransport_,
               VmTransport_, 
               nEcgPoints);
   
   ArrayView<double> ecgs = ecgsTransport_;
   double* ecgsSendBuf=&ecgs[0];
   double ecgsRecvBuf[nEcgPoints];
   MPI_Allreduce(ecgsSendBuf, ecgsRecvBuf, nEcgPoints, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);  

   for(int ii=0; ii<ecgs.size(); ii++)
   {
	ecgs[ii]=ecgsRecvBuf[ii];
   }
 
   stopTimer(sensorEvalTimer);
}

