#include "reactionFactory.hh"

#include <cassert>
#include <cstdlib>
#include <dlfcn.h>
#include <unordered_map>
#include <mpi.h>
#include "object_cc.hh"
#include "ioUtils.h"
#include "mpiUtils.h"

#include "TT04_CellML_Reaction.hh" // TT04 implementation from CellML (Nov 2011)
#include "TT06_CellML_Reaction.hh" // TT06 implementation from CellML (Nov 2011)
#include "TT06Dev_Reaction.hh"
#include "pade.hh"
#include "TT06_RRG_Reaction.hh"    // TT06 with modifications from Rice et al.
#include "OHaraRudy_Reaction.hh"    // OHara Rudy .
//#include "Grandi_Reaction.hh"      // Grandi
#include "ReactionFHN.hh"
#include "NullReaction.hh"
#include "TestReaction.hh"
#include "ConstantReaction.hh"
#include "ThreadServer.hh"
#include "Anatomy.hh"
#include "string.h"

#include <iostream>
using namespace std;

static unordered_map<string,reactionFactoryFunction> g_factoryFromMethodName;

Reaction* reactionFactory(const string& name, double dt, const int numPoints,
                          const ThreadTeam& group)
{
   static bool first = true;
   if (first)
   {
      registerBuiltinReactions();
      first = false;
   }
   
   OBJECT* obj = objectFind(name, "REACTION");
   string method; objectGet(obj, "method", method, "undefined");

   auto iter = g_factoryFromMethodName.find(method);
   if (iter != g_factoryFromMethodName.end())
   {
      return iter->second(obj, dt, numPoints, group);
   }
   string filename = method;
   if (filename[0]!='/')
   {
      filename = "./"+filename;
   }
   if (filetest(filename.c_str(),S_IFREG) == 0)
   {
      //try to load in the factory method
      
      void* handle = dlopen(filename.c_str(), RTLD_NOW|RTLD_LOCAL);
      if (handle)
      {
         Reaction* (*factoryMethod)(OBJECT*,const double,const int,const ThreadTeam&) = reinterpret_cast<Reaction*(*)(OBJECT*,const double,const int,const ThreadTeam&)>(dlsym(handle,"factory"));
         if (factoryMethod)
         {
            return factoryMethod(obj, dt, numPoints, group);
         }
      }
      else
      {
         cerr << "Cant load dynamic module " << filename << ": " << dlerror() << endl;
      }
   }
   {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   
      if (myRank == 0)
         cerr<<"ERROR: Undefined reaction model in reactionFactory"<<endl;
      assert(false); // reachable only due to bad input
      return 0;
   }
}

void registerReactionFactory(const string method, reactionFactoryFunction scanFunc)
{
   g_factoryFromMethodName[method] = scanFunc;
}

