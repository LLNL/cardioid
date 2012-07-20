#ifndef STIMULUS_HH
#define STIMULUS_HH

#include "VectorDouble32.hh"

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
             VectorDouble32& dVmDiffusion)
   {
      if (time > t0_ && time < tf_)
         subClassStim(time, dVmDiffusion);
   }

   virtual void subClassStim(double time,
             VectorDouble32& dVmDiffusion) = 0;

 private:
   double t0_;
   double tf_;
};

#endif
