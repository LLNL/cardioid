#ifndef RANDOM_PULSE_HH
#define RANDOM_PULSE_HH

#include "Pulse.hh"

class RandomPulse:
   public Pulse
{
 public:
   RandomPulse(const double minperiod, const double maxperiod, 
               const double duration,
               const double amplitude, const double phase,
               const double begintime=0.,
               const unsigned short rand_state[3]=0);

   double eval(double time);

 private:
   const double minperiod_;
   const double maxperiod_;
   const double duration_;
   const double amplitude_;
   const double phase_;
   
   double begintime_;
   double endtime_;
   unsigned short rand_state_[3];
};

#endif
