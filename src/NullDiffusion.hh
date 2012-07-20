#ifndef NULL_DIFFUSION_HH
#define NULL_DIFFUSION_HH

#include "Diffusion.hh"

class NullDiffusion : public Diffusion
{
 public:
   NullDiffusion(){};
   
   void updateLocalVoltage(const double* VmLocal){};
   void updateRemoteVoltage(const double* VmRemote){};
   void calc(VectorDouble32& dVm){dVm.assign(dVm.size(), 0.0);};
};


#endif
