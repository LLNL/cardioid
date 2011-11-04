#include "UniformConductivity.hh"

UniformConductivity::UniformConductivity(const UniformConductivityParms& p)
:sigma_(p.sigma)
{
}

   
void UniformConductivity::compute(
   int theta, int phi, SigmaTensorMatrix& sigma)
{
   sigma = sigma_;
}

