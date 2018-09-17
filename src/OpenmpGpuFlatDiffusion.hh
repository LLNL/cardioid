#ifndef OPENMPGPUFLATDIFFUSION_HH
#define OPENMPGPUFLATDIFFUSION_HH

#include "Diffusion.hh"
#include "Anatomy.hh"
#include "LocalGrid.hh"
#include <vector>

class OpenmpGpuFlatDiffusion : public Diffusion
{
 public:
   OpenmpGpuFlatDiffusion(const Anatomy& anatomy, int simLoopType);

   void updateLocalVoltage(ro_larray_ptr<double> VmLocal);
   void updateRemoteVoltage(ro_larray_ptr<double> VmRemote);
   /** omp loop must assign dVm, parallel loop need to increment dVm */
   void calc(rw_larray_ptr<double> dVm);

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
   lazy_array<double> VmBlock_;
   lazy_array<double> dVmBlock_;
   lazy_array<double> sigmaFaceNormal_;
   lazy_array<int> blockFromCell_;
   
   friend void actualCalc(OpenmpGpuFlatDiffusion& self, rw_larray_ptr<double> dVm);
};

void actualCalc(OpenmpGpuFlatDiffusion& self, rw_larray_ptr<double> dVm);
#endif
