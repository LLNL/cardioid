#ifndef STIMULUS_FACTORY
#define STIMULUS_FACTORY

#include <string>
class Stimulus;
class Anatomy;

Stimulus* stimulusFactory(const std::string& name, const Anatomy& anatomy);

#endif
