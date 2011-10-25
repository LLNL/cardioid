#ifndef REACTION_FACTORY_HH
#define REACTION_FACTORY_HH

#include <string>
class Reaction;
class Anatomy;

Reaction* reactionFactory(const std::string& name, const Anatomy& anatomy);

#endif
