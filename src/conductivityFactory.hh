#ifndef CONDUCTIVITY_FACTORY_HH
#define CONDUCTIVITY_FACTORY_HH

#include <string>
class Conductivity;
class Anatomy;

Conductivity* conductivityFactory(const std::string& name, const Anatomy& anatomy);

#endif
