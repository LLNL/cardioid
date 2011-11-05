#ifndef STIMULUS_HH
#define STIMULUS_HH

#include <vector>

class Stimulus
{
 public:
   virtual void stim(double time,
		     std::vector<double>& dVmDiffusion,
		     std::vector<double>& dVmExternal) = 0;
   virtual ~Stimulus(){};
};

#endif
