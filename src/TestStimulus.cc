#include "TestStimulus.hh"
#include <mpi.h>
#include <cmath>

using namespace std;

TestStimulus::TestStimulus(const TestStimulusParms& p)
: targetRank_(p.rank),
  targetCell_(p.cell),
  dVmStim_(-p.iStim),
  tStart_(p.tStart),
  tEnd_(p.tEnd),
  freq_(p.freq)
{
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank_);
}

void TestStimulus::stim(double time,
                        vector<double>& dVmDiffusion,
                        vector<double>& dVmExternal)
{
   if ( myRank_ != targetRank_ )
      return;

   dVmDiffusion[targetCell_] = 0;
   double t = time - floor(time/freq_)*freq_;
   if (t >= tStart_ && t < tEnd_)
      dVmExternal[targetCell_] = dVmStim_;
   else
      dVmExternal[targetCell_] = 0;
}
