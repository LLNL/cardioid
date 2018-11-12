#include "PeriodicPulse.hh"
#include <cmath>
#include <cassert>

using namespace std;


PeriodicPulse::PeriodicPulse(double period, double duration,
                             double amplitude, double phase)
: period_(period),
  duration_(duration),
  amplitude_(amplitude),
  phase_(phase)
{
   // lazy way of checking for bad input.
   assert(period_ > 0);
   assert(phase_ < period_);
   assert(duration_ > 0);
   assert(phase_+duration_ <= period_);
}

double PeriodicPulse::eval(double time) 
{
   double t = time - floor(time/period_)*period_;
   if (t >= phase_ && t < phase_+duration_)
      return amplitude_;
   return 0.0;
}
