#include "BoxStimulus.hh"
#include "Anatomy.hh"
#include <cmath>

using namespace std;

BoxStimulus::BoxStimulus(const BoxStimulusParms& p, const Anatomy& anatomy, Pulse* pulse)
: Stimulus(p.baseParms),
  pulse_(pulse)
{
   // Loop through local points and store the indices of any that are in
   // the stimulated region.
   for (unsigned ii=0; ii<anatomy.nLocal(); ++ii)
   {
      Tuple gt = anatomy.globalTuple(ii);
      if (gt.x() > p.xMin && gt.x() < p.xMax &&
          gt.y() > p.yMin && gt.y() < p.yMax &&
          gt.z() > p.zMin && gt.z() < p.zMax )
         stimList_.push_back(ii);
   }
}

void BoxStimulus::subClassStim(double time,
                               VectorDouble32& dVmDiffusion)
{
   double value = pulse_->eval(time);
   if (value != 0)
      for (unsigned ii=0; ii<stimList_.size(); ++ii)
         dVmDiffusion[stimList_[ii]] += value;
}
