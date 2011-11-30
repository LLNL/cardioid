#include "UniformConductivity.hh"

UniformConductivity::UniformConductivity(const UniformConductivityParms& p)
:sigma_(p.sigma)
{
}

   
void UniformConductivity::compute(const AnatomyCell&, SigmaTensorMatrix& sigma)
{
   sigma = sigma_;
}

SigmaTensorMatrix UniformConductivity::defaultValue()
{
   return sigma_;
}
