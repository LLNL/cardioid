#include "PointStimulus.hh"
#include "Anatomy.hh"

using namespace std;

PointStimulus::PointStimulus(const PointStimulusParms& p, const Anatomy& anatomy)
: Stimulus(p.baseParms),
  pulse_(p.period, p.duration, -p.vStim, p.tStart)
{

  // loop through grid points on this task, check against target cell
  cellLocal_ = false;
  for (unsigned jj=0; jj<anatomy.size(); ++jj)
  {
    if (anatomy.gid(jj) == p.cell)
    {        
      cellLocal_ = true;
      localInd_ = jj;
      break;
    }
  }
}

void PointStimulus::subClassStim(double time,
                                 vector<double>& dVmDiffusion,
                                 vector<double>& dVmExternal)
{
   if (cellLocal_)
   {
      dVmExternal[localInd_] += pulse_.eval(time);
   }
}
