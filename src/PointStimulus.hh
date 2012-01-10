#ifndef POINT_STIMULUS_HH
#define POINT_STIMULUS_HH

#include "Stimulus.hh"
#include "PeriodicPulse.hh"

class Anatomy;

struct PointStimulusParms
{
   int cell;
   double vStim;
   double tStart;
   double period;
   double duration;
   StimulusBaseParms baseParms;   
};

class PointStimulus : public Stimulus
{
 public:
   PointStimulus(const PointStimulusParms& p, const Anatomy& anatomy);
   void subClassStim(double time,
                     std::vector<double>& dVmDiffusion,
                     std::vector<double>& dVmExternal);
   
 private:
   int targetCell_;
   bool cellLocal_;
   int localInd_;
   PeriodicPulse pulse_;
};

#endif
