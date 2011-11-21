#ifndef POINT_STIMULUS_HH
#define POINT_STIMULUS_HH

#include "Stimulus.hh"

class Anatomy;

struct PointStimulusParms
{
   int cell;
   double iStim;
   double tStart;
   double tEnd;
   double freq;
};

class PointStimulus : public Stimulus
{
 public:
   PointStimulus(const PointStimulusParms& p, const Anatomy& anatomy);
   void stim(double time,
	     std::vector<double>& dVmDiffusion,
	     std::vector<double>& dVmExternal);
   
 private:
   int targetCell_;
   double dVmStim_;
   double tStart_;
   double tEnd_;
   double freq_;
   bool cellLocal_;
   int localInd_;
};

#endif
