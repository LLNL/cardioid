#ifndef REACTION_FACTORY_HH
#define REACTION_FACTORY_HH
#include <vector>
#include <string>
class ThreadTeam;
class Reaction;
class Anatomy;
using std::string;
using std::vector;

Reaction* reactionFactory(const std::string& name, double dt, const int numPoints,
                          const ThreadTeam &group, const vector<string>& scaleCurrents);

#endif
