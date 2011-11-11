#include "reactionFactory.hh"

#include <cassert>
#include "object_cc.hh"

#include "TT04_bbReaction.hh"         // TT04 implementation from BlueBeats
#include "TT04_CellML_Reaction.hh"    // TT04 implementation from CellML (Nov 2011)
#include "TT06_CellML_Reaction.hh"    // TT06 implementation from CellML (Nov 2011)
using namespace std;

namespace
{
   Reaction* scanTT04_bb(OBJECT* obj, const Anatomy& anatomy);
   Reaction* scanTT04_CellML(OBJECT* obj, const Anatomy& anatomy);
   Reaction* scanTT06_CellML(OBJECT* obj, const Anatomy& anatomy);
}


Reaction* reactionFactory(const string& name, const Anatomy& anatomy)
{
   OBJECT* obj = objectFind(name, "REACTION");
   string method; objectGet(obj, "method", method, "undefined");

   if (method == "undefined")
      assert(1==0);
   else if (method == "TT04_bb" || method == "tenTusscher04_bb")
      return scanTT04_bb(obj, anatomy);
   else if (method == "TT04_CellML" || method == "tenTusscher04_CellML")
      return scanTT04_CellML(obj, anatomy);
   else if (method == "TT06_CellML" || method == "tenTusscher06_CellML")
      return scanTT06_CellML(obj, anatomy);
   
   assert(false); // reachable only due to bad input
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
   Reaction* scanTT04_CellML(OBJECT* obj, const Anatomy& anatomy)
   {
      TT04_CellML_Reaction::IntegratorType integrator;
      string tmp;
      objectGet(obj, "integrator", tmp, "rushLarson");
      if      (tmp == "rushLarson")   integrator = TT04_CellML_Reaction::rushLarson;
      else if (tmp == "forwardEuler") integrator = TT04_CellML_Reaction::forwardEuler;
      else    assert(false);    
      
      return new TT04_CellML_Reaction(anatomy, integrator);
   }
}

namespace
{
   Reaction* scanTT06_CellML(OBJECT* obj, const Anatomy& anatomy)
   {
      return new TT06_CellML_Reaction(anatomy);
   }
}
