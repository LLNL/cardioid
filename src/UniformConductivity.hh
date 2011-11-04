#ifndef UNIFORM_CONDUCTIVTY_HH
#define UNIFORM_CONDUCTIVTY_HH

#include "Conductivity.hh"

struct UniformConductivityParms
{
   SigmaTensorMatrix sigma;
};


class UniformConductivity : public Conductivity
{
 public:
   UniformConductivity(const UniformConductivityParms& p);
   
   void compute(int theta, int phi, SigmaTensorMatrix& sigma);
 private:
   SigmaTensorMatrix sigma_;
};

#endif
