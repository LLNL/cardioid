#include "PointStimulus.hh"
#include "Anatomy.hh"

using namespace std;

PointStimulus::PointStimulus(const PointStimulusParms& p,
                             const Anatomy& anatomy,
                             Pulse* pulse)
: Stimulus(p.baseParms),
  pulse_(pulse)
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

int PointStimulus::subClassStim(double time,
                                rw_mgarray_ptr<double> _dVmDiffusion)
{
   rw_array_ptr<double> dVmDiffusion = _dVmDiffusion.useOn(CPU);
   if (cellLocal_)
   {
      dVmDiffusion[localInd_] += pulse_->eval(time);
      return 1;
   }
   return 0;
}

int PointStimulus::nStim()
{
   if (cellLocal_) 
      return 1;
   return 0;
}
