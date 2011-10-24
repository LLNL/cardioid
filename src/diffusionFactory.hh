#ifndef DIFFUSION_FACTORY_HH
#define DIFFUSION_FACTORY_HH

#include <string>
class Diffusion;
class Anatomy;

Diffusion* diffusionFactory(const std::string& name, const Anatomy& anatomy);

#endif
