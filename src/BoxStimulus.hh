#ifndef BOX_STIMULUS_HH
#define BOX_STIMULUS_HH

#include "Stimulus.hh"
#include "Pulse.hh"
#include "TransportCoordinator.hh"
#include <string>

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
               Pulse* pulse, const std::string& name);
   int subClassStim(double time,
		    Managed<ArrayView<double>> dVmDiffusion);
   int nStim();
   
 private:
   Pulse* pulse_;
   TransportCoordinator<PinnedVector<int> > stimListTransport_;
};

#endif
