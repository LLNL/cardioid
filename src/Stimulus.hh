#ifndef STIMULUS_HH
#define STIMULUS_HH

#include <vector>

struct StimulusBaseParms
{
   double t0;
   double tf;
};


class Stimulus
{
 public:
   Stimulus(const StimulusBaseParms& p)
   : t0_(p.t0), tf_(p.tf){}
   
   virtual ~Stimulus(){};

   void stim(double time,
             std::vector<double>& dVmDiffusion,
             std::vector<double>& dVmExternal)
   {
      if (time > t0_ && time < tf_)
         subClassStim(time, dVmDiffusion, dVmExternal);
   }

   virtual void subClassStim(double time,
             std::vector<double>& dVmDiffusion,
             std::vector<double>& dVmExternal) = 0;

 private:
   double t0_;
   double tf_;
};

#endif
