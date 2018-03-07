#ifndef TEST_STIMULUS_HH
#define TEST_STIMULUS_HH

#include "Stimulus.hh"
#include "Pulse.hh"

struct TestStimulusParms
{
   int rank;
   int cell;
   StimulusBaseParms baseParms;
};

class TestStimulus : public Stimulus
{
 public:
   TestStimulus(const TestStimulusParms& p, Pulse* pulse);
   int subClassStim(double time,
		    Managed<ArrayView<double>> dVmDiffusion);
   int nStim();
   
 private:
   int targetRank_;
   int targetCell_;
   int myRank_;
   Pulse* pulse_;
};

#endif
