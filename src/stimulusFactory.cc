#include "stimulusFactory.hh"

#include <cassert>
#include "object_cc.hh"
#include "Stimulus.hh"
#include "PointStimulus.hh"
#include "TestStimulus.hh"

using namespace std;

namespace
{
   Stimulus* scanPointStimulus(OBJECT* obj, const Anatomy& anatomy);
   Stimulus* scanTestStimulus(OBJECT* obj);
}


Stimulus* stimulusFactory(const std::string& name, const Anatomy& anatomy)
{
   OBJECT* obj = objectFind(name, "STIMULUS");
   string method; objectGet(obj, "method", method, "undefined");

   if (method == "undefined")
      assert(false);
   else if (method == "point")
     return scanPointStimulus(obj, anatomy);
   else if (method == "test")
      return scanTestStimulus(obj);

   assert(false); // reachable only due to bad input
}


namespace
{
   Stimulus* scanPointStimulus(OBJECT* obj, const Anatomy& anatomy)
   {
      PointStimulusParms p;
      objectGet(obj, "cell",   p.cell,   "0");
      objectGet(obj, "freq",   p.freq,   "1000");
      objectGet(obj, "iStim",  p.iStim,  "-52");
      objectGet(obj, "tEnd",   p.tEnd,   "1");
      objectGet(obj, "tStart", p.tStart, "2");
      return new PointStimulus(p, anatomy);
   }
   Stimulus* scanTestStimulus(OBJECT* obj)
   {
      TestStimulusParms p;
      objectGet(obj, "cell",   p.cell,   "0");
      objectGet(obj, "freq",   p.freq,   "1000");
      objectGet(obj, "iStim",  p.iStim,  "-52");
      objectGet(obj, "rank",   p.rank,   "0");
      objectGet(obj, "tEnd",   p.tEnd,   "1");
      objectGet(obj, "tStart", p.tStart, "2");
      return new TestStimulus(p);
   }
}
