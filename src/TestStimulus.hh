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
   void subClassStim(double time,
                     std::vector<double>& dVmDiffusion);
   
 private:
   int targetRank_;
   int targetCell_;
   int myRank_;
   Pulse* pulse_;
};

#endif
