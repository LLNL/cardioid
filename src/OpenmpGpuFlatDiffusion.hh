#ifndef OPENMPGPUFLATDIFFUSION_HH
#define OPENMPGPUFLATDIFFUSION_HH

#include "Diffusion.hh"
#include "TransportCoordinator.hh"
#include "Anatomy.hh"
#include "LocalGrid.hh"
#include <vector>

class OpenmpGpuFlatDiffusion : public Diffusion
{
 public:
   OpenmpGpuFlatDiffusion(const Anatomy& anatomy, int simLoopType);

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
   int nx_;
   int ny_;
   int nz_;
   bool skip_reorder_;
   
   double Vm_;
   double dVm_;
   LocalGrid localGrid_;
   TransportCoordinator<PinnedVector<double> > VmBlock_;
   TransportCoordinator<PinnedVector<double> > dVmBlock_;
   TransportCoordinator<PinnedVector<double> > sigmaFaceNormal_;
   TransportCoordinator<PinnedVector<int> > blockFromCell_;
   
   friend void actualCalc(OpenmpGpuFlatDiffusion& self, ArrayView<double> dVm);
};

void actualCalc(OpenmpGpuFlatDiffusion& self, ArrayView<double> dVm);
#endif
