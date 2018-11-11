#ifndef POINT_STIMULUS_HH
#define POINT_STIMULUS_HH

#include "Stimulus.hh"
#include "Long64.hh"
#include "Pulse.hh"

class Anatomy;

struct PointStimulusParms
{
   Long64 cell;
   StimulusBaseParms baseParms;   
};

class PointStimulus : public Stimulus
{
 public:
   PointStimulus(const PointStimulusParms& p, const Anatomy& anatomy, 
                 Pulse* pulse);
   int subClassStim(double time,
		    rw_mgarray_ptr<double> dVmDiffusion);
   int nStim();
   
 private:
   Long64 targetCell_;
   bool cellLocal_;
   int localInd_;
   Pulse* pulse_;
};

#endif
