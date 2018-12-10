#ifndef BOX_STIMULUS_HH
#define BOX_STIMULUS_HH

#include "Stimulus.hh"
#include "Pulse.hh"
#include "lazy_array.hh"
#include <string>

class Anatomy;

struct BoxStimulusParms
{
   double period;
   bool isXYZ;
   double x;
   double y;
   double z;
   double length;
   // TO DO: ELIMINATE GID RANGES
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
		    rw_mgarray_ptr<double> dVmDiffusion);
   int nStim();

 private:
   Pulse* pulse_;
   lazy_array<int> stimListTransport_;
};

#endif
