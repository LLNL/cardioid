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

   void updateLocalVoltage(const Managed<ArrayView<double>> VmLocal);
   void updateRemoteVoltage(const Managed<ArrayView<double>> VmRemote);
   /** omp loop must assign dVm, parallel loop need to increment dVm */
   void calc(Managed<ArrayView<double>> dVm);

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
   TransportCoordinator<PinnedVector<double> > VmBlock_;
   TransportCoordinator<PinnedVector<double> > sigmaFaceNormal_;
   TransportCoordinator<PinnedVector<int> > blockFromRed_;
   TransportCoordinator<PinnedVector<int> > cellFromRed_;
   TransportCoordinator<PinnedVector<int> > cellLookup_;
   
   friend void actualCalc(OpenmpGpuRedblackDiffusion& self, ArrayView<double> dVm);
};

void actualCalc(OpenmpGpuRedblackDiffusion& self, ArrayView<double> dVm);
#endif
