#ifndef REACTION_FACTORY_HH
#define REACTION_FACTORY_HH

#include <string>
class ThreadTeam;
class Reaction;
class Anatomy;

Reaction* reactionFactory(const std::string& name, double dt, Anatomy& anatomy,
                          const ThreadTeam &group);

#endif
