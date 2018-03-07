#ifndef NULL_DIFFUSION_HH
#define NULL_DIFFUSION_HH

#include "Diffusion.hh"
#include <vector>
#include "Anatomy.hh"

class NullDiffusion : public Diffusion
{
 public:
   NullDiffusion(const Anatomy& anatomy, int simLoopType)
   : simLoopType_(simLoopType),
     blockIndex_(anatomy.size(), 0),
     Vm_(0),
     dVm_(0)
   {};
   
   void updateLocalVoltage(const Managed<ArrayView<double>> VmLocal) {}
   void updateRemoteVoltage(const Managed<ArrayView<double>> VmRemote) {}
   /** omp loop must assign dVm, parallel loop need to increment dVm */
   void calc(Managed<ArrayView<double>> dVm)
   {
      if (simLoopType_ == 0)
      {
         ArrayView<double> dVmHost = dVm.modifyOnHost();
         for (int ii=0; ii<dVmHost.size(); ++ii)
         {
            dVmHost[ii] = 0;
         }
      }
   };

   unsigned* blockIndex() {return &blockIndex_[0];}
   double* VmBlock() {return &Vm_;}
   double* dVmBlock() {return &dVm_;}
   
 private:
   int simLoopType_;
   std::vector<unsigned> blockIndex_;
   double Vm_;
   double dVm_;
   
};


#endif
