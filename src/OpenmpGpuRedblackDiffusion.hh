#ifndef OPENMPGPUREDBLACKDIFFUSION_HH
#define OPENMPGPUREDBLACKDIFFUSION_HH

#include "Diffusion.hh"
#include "TransportCoordinator.hh"
#include "Anatomy.hh"
#include "LocalGrid.hh"
#include <vector>

class OpenmpGpuRedblackDiffusion : public Diffusion
{
 public:
   OpenmpGpuRedblackDiffusion(const Anatomy& anatomy, int simLoopType);

   void updateLocalVoltage(const double* VmLocal);
   void updateRemoteVoltage(const double* VmRemote);
   /** omp loop must assign dVm, parallel loop need to increment dVm */
   void calc(VectorDouble32& dVm);

   unsigned* blockIndex();
   double* VmBlock();
   double* dVmBlock();
   
 private:
   int simLoopType_;
   int nLocal_;
   int nRemote_;
   int nCells_;
   int nRed_;
   int nBlack_;
   int nRedLocal_;
   int nx_;
   int ny_;
   int nz_;
   
   double Vm_;
   double dVm_;
   LocalGrid localGrid_;
   TransportCoordinator<std::vector<double> > VmBlock_;
   TransportCoordinator<std::vector<double> > sigmaFaceNormal_;
   TransportCoordinator<std::vector<int> > blockFromRed_;
   TransportCoordinator<std::vector<int> > cellFromRed_;
   TransportCoordinator<std::vector<int> > cellLookup_;
   
   friend void actualCalc(OpenmpGpuRedblackDiffusion& self, VectorDouble32& dVm);
};

void actualCalc(OpenmpGpuRedblackDiffusion& self, VectorDouble32& dVm);
#endif
