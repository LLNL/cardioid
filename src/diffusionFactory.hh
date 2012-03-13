#ifndef DIFFUSION_FACTORY_HH
#define DIFFUSION_FACTORY_HH

#include <string>
class Diffusion;
class Anatomy;
class CoreGroup;

Diffusion* diffusionFactory(const std::string& name, const Anatomy& anatomy, const CoreGroup& threadInfo);

#endif
