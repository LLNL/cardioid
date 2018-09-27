#ifndef CUDA_DIFFUSION_HH
#define CUDA_DIFFUSION_HH

#include "Diffusion.hh"
#include "lazy_array.hh"
#include "Anatomy.hh"
#include "LocalGrid.hh"
#include <vector>

class CUDADiffusion : public Diffusion
{
 public:
   CUDADiffusion(const Anatomy& anatomy, int simLoopType, double diffusionScale);

   void updateLocalVoltage(ro_mgarray_ptr<double> VmLocal);
   void updateRemoteVoltage(ro_mgarray_ptr<double> VmRemote);
   /** omp loop must assign dVm, parallel loop need to increment dVm */
   void calc(rw_mgarray_ptr<double> dVm);

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
   lazy_array<double> VmBlock_;
   lazy_array<double> dVmBlock_;
   lazy_array<double> sigmaFaceNormal_;
   lazy_array<int> cellLookup_;
   
};

#endif
