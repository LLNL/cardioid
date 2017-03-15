#ifndef GPU_DIFFUSION_HH
#define GPU_DIFFUSION_HH

#include "Diffusion.hh"
#include <vector>
#include "Anatomy.hh"

class GPUDiffusion : public Diffusion
{
 public:
   GPUDiffusion(const Anatomy& anatomy, int simLoopType);

   void updateLocalVoltage(const double* VmLocal);
   void updateRemoteVoltage(const double* VmRemote);
   /** omp loop must assign dVm, parallel loop need to increment dVm */
   void calc(VectorDouble32& dVm);

   unsigned* blockIndex();
   double* VmBlock();
   double* dVmBlock();
   
 private:
   int simLoopType_;
   std::vector<unsigned> blockIndex_;
   double Vm_;
   double dVm_;
   
};


#endif
