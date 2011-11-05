#ifndef TEST_STIMULUS_HH
#define TEST_STIMULUS_HH

#include "Stimulus.hh"

struct TestStimulusParms
{
   int rank;
   int cell;
   double iStim;
   double tStart;
   double tEnd;
   double freq;
};

class TestStimulus : public Stimulus
{
 public:
   TestStimulus(const TestStimulusParms& p);
   void stim(double time,
	     std::vector<double>& dVmDiffusion,
	     std::vector<double>& dVmExternal);
   
 private:
   int targetRank_;
   int targetCell_;
   int myRank_;
   double dVmStim_;
   double tStart_;
   double tEnd_;
   double freq_;
};

#endif
