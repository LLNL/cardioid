#include "reactionFactory.hh"

#include <cassert>
#include "object_cc.hh"

#include "TT04_bbReaction.hh" // TT04 implementation from BlueBeats
#include "TT04Reaction.hh"    // 
			      // new TT04 implementation from Matthias (Oct 2011)
using namespace std;

namespace
{
   Reaction* scanTT04_bb(OBJECT* obj, const Anatomy& anatomy);
   Reaction* scanTT04(OBJECT* obj, const Anatomy& anatomy);
}


Reaction* reactionFactory(const string& name, const Anatomy& anatomy)
{
   OBJECT* obj = objectFind(name, "REACTION");
   string method; objectGet(obj, "method", method, "undefined");

   if (method == "undefined")
      assert(1==0);
   else if (method == "TT04_bb" || method == "tenTusscher04_bb")
      return scanTT04_bb(obj, anatomy);
   else if (method == "TT04" || method == "tenTusscher04")
      return scanTT04(obj, anatomy);
   
   assert(1==0); // reachable only due to bad input
}

namespace
{
   Reaction* scanTT04_bb(OBJECT* obj, const Anatomy& anatomy)
   {
      return new TT04_bbReaction(anatomy);
   }
}

namespace
{
   Reaction* scanTT04(OBJECT* obj, const Anatomy& anatomy)
   {
      return new TT04Reaction(anatomy);
   }
}
