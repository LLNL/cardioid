#ifndef OPENMPGPUREDBLACKDIFFUSION_HH
#define OPENMPGPUREDBLACKDIFFUSION_HH

#include "Diffusion.hh"
#include "Anatomy.hh"
#include "LocalGrid.hh"
#include <vector>

class OpenmpGpuRedblackDiffusion : public Diffusion
{
 public:
   OpenmpGpuRedblackDiffusion(const Anatomy& anatomy, int simLoopType);

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
   int nRed_;
   int nBlack_;
   int nRedLocal_;
   int nx_;
   int ny_;
   int nz_;
   
   //double Vm_;
   //double dVm_;
   LocalGrid localGrid_;
   lazy_array<double> VmBlock_;
   lazy_array<double> sigmaFaceNormal_;
   lazy_array<int> blockFromRed_;
   lazy_array<int> cellFromRed_;
   lazy_array<int> cellLookup_;
   
   friend void actualCalc(OpenmpGpuRedblackDiffusion& self, rw_larray_ptr<double> dVm);
};

void actualCalc(OpenmpGpuRedblackDiffusion& self, rw_larray_ptr<double> dVm);
#endif
