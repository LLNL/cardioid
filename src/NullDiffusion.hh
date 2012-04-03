#ifndef NULL_DIFFUSION_HH
#define NULL_DIFFUSION_HH

#include "Diffusion.hh"

class NullDiffusion : public Diffusion
{
 public:
   NullDiffusion(){};
   
   void updateLocalVoltage(const double* VmLocal){};
   void updateRemoteVoltage(const double* VmRemote){};
   void calc(std::vector<double>& dVm){};
};


#endif
