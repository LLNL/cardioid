#ifndef PERIODIC_PULSE_HH
#define PERIODIC_PULSE_HH

#include "Pulse.hh"

/** Periodic square pulse useful for creating external stimulus
 * functions.
 *
 *  let t = time - floor(time/period)*period.
 *  when t >= phase && t < phase_duration eval returns amplitude.
 *  returns zero otherwise.
*/
class PeriodicPulse:
   public Pulse
{
 public:
   PeriodicPulse(double period, double duration,
                 double amplitude, double phase);

   double eval(double time);

 private:
   double period_;
   double duration_;
   double amplitude_;
   double phase_;
};

#endif
