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
   
   void updateLocalVoltage(ro_larray_ptr<double> VmLocal) {};
   void updateRemoteVoltage(ro_larray_ptr<double> VmRemote) {};
   /** omp loop must assign dVm, parallel loop need to increment dVm */
   void calc(rw_larray_ptr<double> dVm)
   {
      if (simLoopType_ == 0)
      {
         ContextRegion region(CPU);
         dVm.use();
         for (int ii=0; ii<dVm.size(); ++ii)
         {
            dVm[ii] = 0;
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
