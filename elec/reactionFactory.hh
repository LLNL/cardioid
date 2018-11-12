#ifndef REACTION_FACTORY_HH
#define REACTION_FACTORY_HH
#include <vector>
#include <string>
#include "object.h"
class ThreadTeam;
class Reaction;

typedef Reaction* (*reactionFactoryFunction)(OBJECT* obj, const double dt, const int numPoints, const ThreadTeam& group);

Reaction* reactionFactory(const std::string& name, double dt, const int numPoints,
                          const ThreadTeam &group);

void registerReactionFactory(const std::string method, reactionFactoryFunction scanFunc);

void registerBuiltinReactions();

#ifdef DYNAMIC_REACTION
#define REACTION_FACTORY(name) extern "C" Reaction* factory
#define FRIEND_FACTORY(name) friend Reaction* ::factory
#else
#define REACTION_FACTORY(name) Reaction* reactionFactoryFor##name
#define FRIEND_FACTORY(name) friend Reaction* ::reactionFactoryFor##name
#endif

#endif
