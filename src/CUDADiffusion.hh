#ifndef CUDA_DIFFUSION_HH
#define CUDA_DIFFUSION_HH

#include "Diffusion.hh"
#include "TransportCoordinator.hh"
#include "Anatomy.hh"
#include "LocalGrid.hh"
#include <vector>

class CUDADiffusion : public Diffusion
{
 public:
   CUDADiffusion(const Anatomy& anatomy, int simLoopType);

   void updateLocalVoltage(const double* VmLocal);
   void updateRemoteVoltage(const double* VmRemote);
   /** omp loop must assign dVm, parallel loop need to increment dVm */
   void calc(VectorDouble32& dVm);

   unsigned* blockIndex();
   double* VmBlock();
   double* dVmBlock();
   
 private:
   const Anatomy& anatomy_; 
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
   TransportCoordinator<std::vector<double> > dVmBlock_;
   TransportCoordinator<std::vector<double> > sigmaFaceNormal_;
   TransportCoordinator<std::vector<int> > cellLookup_;
   
   friend void actualCalc(CUDADiffusion& self, VectorDouble32& dVm);
};

void actualCalc(CUDADiffusion& self, VectorDouble32& dVm);
#endif
