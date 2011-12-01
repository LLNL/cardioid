#include "BoxStimulus.hh"
#include "Anatomy.hh"
#include <cmath>

using namespace std;

BoxStimulus::BoxStimulus(const BoxStimulusParms& p, const Anatomy& anatomy)
: dVmStim_(-p.iStim),
  tStart_(p.tStart),
  freq_(p.freq),
  duration_(p.duration)
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

void BoxStimulus::stim(double time,
		       vector<double>& dVmDiffusion,
		       vector<double>& dVmExternal)
{
   double trel = time - tStart_;  // stimulus starts at tStart_
   double t = trel - floor(trel/freq_)*freq_;  
   if (t < duration_)
      for (unsigned ii=0; ii<stimList_.size(); ++ii)
	 dVmExternal[stimList_[ii]] = dVmStim_;
}
