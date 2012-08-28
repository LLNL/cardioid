#ifndef NULL_DIFFUSION_HH
#define NULL_DIFFUSION_HH

#include "Diffusion.hh"
#include <vector>
#include "Anatomy.hh"

class NullDiffusion : public Diffusion
{
 public:
   NullDiffusion(const Anatomy& anatomy)
   : blockIndex_(anatomy.size(), 0),
     Vm_(0),
     dVm_(0)
   {};
   
   void updateLocalVoltage(const double* VmLocal){};
   void updateRemoteVoltage(const double* VmRemote){};
   void calc(VectorDouble32& dVm){dVm.assign(dVm.size(), 0.0);};

   unsigned* blockIndex() {return &blockIndex_[0];}
   double* VmBlock() {return &Vm_;}
   double* dVmBlock() {return &dVm_;}
   
 private:
   std::vector<unsigned> blockIndex_;
   double Vm_;
   double dVm_;
   
};


#endif
