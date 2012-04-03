#ifndef DIFFUSION_HH
#define DIFFUSION_HH

#include <vector>

class Diffusion
{
 public:
   virtual ~Diffusion(){};
   virtual void
   calc(const std::vector<double>& Vm, std::vector<double>& dVm, double *recv_buf, int nLocal) = 0;
};

#endif
