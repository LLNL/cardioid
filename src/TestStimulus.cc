#include "TestStimulus.hh"
#include <mpi.h>
#include <cmath>

using namespace std;

TestStimulus::TestStimulus(const TestStimulusParms& p, Pulse* pulse)
: Stimulus(p.baseParms),
  targetRank_(p.rank),
  targetCell_(p.cell),
  pulse_(pulse)
{
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank_);
}

int TestStimulus::subClassStim(double time,
                                VectorDouble32& dVmDiffusion)
{
   if ( myRank_ != targetRank_ )
      return 0;
   
   dVmDiffusion[targetCell_] = pulse_->eval(time);
   return 1;
}

int TestStimulus::nStim()
{
   if (myRank_ == targetRank_)
      return 1;
   return 0;
}
