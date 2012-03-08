#ifndef NULL_DIFFUSION_HH
#define NULL_DIFFUSION_HH

#include "Diffusion.hh"

class NullDiffusion : public Diffusion
{
 public:
   NullDiffusion(){};
   
   void calc(const std::vector<double>& Vm, std::vector<double>& dVm, double *recv_buf, int nLocal){};
   void calc_simd(const std::vector<double>& Vm, std::vector<double>& dVm, double *recv_buf, int nLocal){};
};


#endif
