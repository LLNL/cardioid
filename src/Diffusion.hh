#ifndef DIFFUSION_HH
#define DIFFUSION_HH

#include <vector>

class Diffusion
{
 public:
   virtual void diffusion(
      const std::vector<double>& Vm,
      std::vector<double>& Istim) = 0;
};


#endif
