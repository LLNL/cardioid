#ifndef BOX_STIMULUS_HH
#define BOX_STIMULUS_HH

#include "Stimulus.hh"

class Anatomy;

struct BoxStimulusParms
{
   double iStim;
   double tStart;
   double freq;
   double duration;
   double xMin;
   double yMin;
   double zMin;
   double xMax;
   double yMax;
   double zMax;
};

class BoxStimulus : public Stimulus
{
 public:
   BoxStimulus(const BoxStimulusParms& p, const Anatomy& anatomy);
   void stim(double time,
             std::vector<double>& dVmDiffusion,
             std::vector<double>& dVmExternal);
   
 private:
   double dVmStim_;
   double tStart_;
   double freq_;
   double duration_;
   std::vector<int> stimList_;
};

#endif
