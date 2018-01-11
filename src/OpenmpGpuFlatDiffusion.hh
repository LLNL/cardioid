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
   int nx_;
   int ny_;
   int nz_;
   bool skip_reorder_;
   
   double Vm_;
   double dVm_;
   LocalGrid localGrid_;
   TransportCoordinator<std::vector<double> > VmBlock_;
   TransportCoordinator<std::vector<double> > dVmBlock_;
   TransportCoordinator<std::vector<double> > sigmaFaceNormal_;
   TransportCoordinator<std::vector<int> > blockFromCell_;
   
   friend void actualCalc(OpenmpGpuFlatDiffusion& self, VectorDouble32& dVm);
};

void actualCalc(OpenmpGpuFlatDiffusion& self, VectorDouble32& dVm);
#endif
