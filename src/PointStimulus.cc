#include "PointStimulus.hh"
#include "Anatomy.hh"
#include <mpi.h>
#include <cmath>
#include <iostream>
using namespace std;

PointStimulus::PointStimulus(const PointStimulusParms& p, const Anatomy& anatomy)
    : dVmStim_(-p.iStim),tStart_(p.tStart),tStop_(p.tStop),freq_(p.freq), duration_(p.duration)
{

  // loop through grid points on this task, check against target cell
  cellLocal_ = false;
  const unsigned gid = p.cell;
  for (unsigned jj=0; jj<anatomy.size(); ++jj)
  {
    if (anatomy.gid(jj) == gid)
    {        
      cellLocal_ = true;
      localInd_ = jj;
    }
  }
}

void PointStimulus::stim(double time,
                        vector<double>& dVmDiffusion,
                        vector<double>& dVmExternal)
{
  // apply stimulus only within tStart_ to tStop_ window
  if (time >= tStart_ && (time <= tStop_ || tStop_ <= 0.0))
  {
    if (cellLocal_)
    {
      //dVmDiffusion[localInd_] = 0.;  //ewd:  should this be zeroed?
      
      double trel = time - tStart_;  // stimulus starts at tStart_

      double t = trel - floor(trel/freq_)*freq_;  
      if (t < duration_)
        dVmExternal[localInd_] = dVmStim_;
      else
        dVmExternal[localInd_] = 0.;
    }
  }
}
