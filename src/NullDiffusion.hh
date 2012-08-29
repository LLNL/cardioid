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
   
   void updateLocalVoltage(const double* VmLocal){};
   void updateRemoteVoltage(const double* VmRemote){};
   /** omp loop must assign dVm, parallel loop need to increment dVm */
   void calc(VectorDouble32& dVm)
   {
      if (simLoopType_ == 0)
         dVm.assign(dVm.size(), 0.0);
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
