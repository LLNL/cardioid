#ifndef CONDUCTIVITY_FACTORY_HH
#define CONDUCTIVITY_FACTORY_HH

#include <string>
class Conductivity;

Conductivity* conductivityFactory(const std::string& name);

#endif
