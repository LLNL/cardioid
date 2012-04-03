#ifndef DIFFUSION_HH
#define DIFFUSION_HH

#include <vector>

class Diffusion
{
 public:
   virtual ~Diffusion(){};
   virtual void updateLocalVoltage(const double* VmLocal) = 0;
   virtual void updateRemoteVoltage(const double* VmRemote) = 0;
   virtual void calc(std::vector<double>& dVm) = 0;
};

#endif
