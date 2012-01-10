#ifndef TEST_STIMULUS_HH
#define TEST_STIMULUS_HH

#include "Stimulus.hh"
#include "PeriodicPulse.hh"

struct TestStimulusParms
{
   int rank;
   int cell;
   double vStim;
   double tStart;
   double duration;
   double period;
   StimulusBaseParms baseParms;
};

class TestStimulus : public Stimulus
{
 public:
   TestStimulus(const TestStimulusParms& p);
   void subClassStim(double time,
                     std::vector<double>& dVmDiffusion,
                     std::vector<double>& dVmExternal);
   
 private:
   int targetRank_;
   int targetCell_;
   int myRank_;
   PeriodicPulse pulse_;
};

#endif
