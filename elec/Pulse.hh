#ifndef PULSE_HH
#define PULSE_HH

class Pulse
{
 public:
   Pulse(){};

   virtual double eval(double time)=0;
};

#endif
