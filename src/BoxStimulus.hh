#ifndef BOX_STIMULUS_HH
#define BOX_STIMULUS_HH

#include "Stimulus.hh"
#include "PeriodicPulse.hh"

class Anatomy;

struct BoxStimulusParms
{
   double vStim;
   double tStart;
   double period;
   double duration;
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
   BoxStimulus(const BoxStimulusParms& p, const Anatomy& anatomy);
   void subClassStim(double time,
                     std::vector<double>& dVmDiffusion);
   
 private:
   PeriodicPulse pulse_;
   std::vector<int> stimList_;
};

#endif
