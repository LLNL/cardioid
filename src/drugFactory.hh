#ifndef DRUG_FACTORY
#define DRUG_FACTORY

#include <string>
class Drug;
class Simulate;

Drug* drugFactory(const std::string& name, const Simulate& sim);

#endif
