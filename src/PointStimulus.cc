#include "PointStimulus.hh"
#include "Anatomy.hh"
#include <mpi.h>
#include <cmath>

using namespace std;

PointStimulus::PointStimulus(const PointStimulusParms& p, const Anatomy& anatomy)
    : dVmStim_(-p.iStim),tStart_(p.tStart),tEnd_(p.tEnd),freq_(p.freq)
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
  if (cellLocal_)
  {
    dVmDiffusion[localInd_] = 0.;
    double t = time - floor(time/freq_)*freq_;
    if (t >= tStart_ && t < tEnd_)
      dVmExternal[localInd_] = dVmStim_;
    else
      dVmExternal[localInd_] = 0.;
  }
}
