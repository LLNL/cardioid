#ifndef BOX_STIMULUS_HH
#define BOX_STIMULUS_HH

#include "Stimulus.hh"
#include "Pulse.hh"

class Anatomy;

struct BoxStimulusParms
{
   double period;
   double xMin;
   double yMin;
   double zMin;
   double xMax;
   double yMax;
   double zMax;
   StimulusBaseParms baseParms;
};

class BoxStimulus : public Stimulus
{
 public:
   BoxStimulus(const BoxStimulusParms& p, const Anatomy& anatomy,
               Pulse* pulse);
   void subClassStim(double time,
                     std::vector<double>& dVmDiffusion);
   
 private:
   Pulse* pulse_;
   std::vector<int> stimList_;
};

#endif
