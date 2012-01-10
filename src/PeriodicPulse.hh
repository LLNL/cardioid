#ifndef PERIODIC_PULSE_HH
#define PERIODIC_PULSE_HH

/** Periodic square pulse useful for creating external stimulus
 * functions.
 *
 *  let t = time - floor(time/period)*period.
 *  when t >= phase && t < phase_duration eval returns amplitude.
 *  returns zero otherwise.
*/
class PeriodicPulse
{
 public:
   PeriodicPulse(double period, double duration,
                 double amplitude, double phase);

   double eval(double time) const;

 private:
   double period_;
   double duration_;
   double amplitude_;
   double phase_;
};

#endif
