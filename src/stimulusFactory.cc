#include "stimulusFactory.hh"

#include <cassert>
#include <sstream>
#include "object_cc.hh"
#include "Anatomy.hh"
#include "Stimulus.hh"
#include "PointStimulus.hh"
#include "TestStimulus.hh"
#include "BoxStimulus.hh"

using namespace std;

namespace
{
   Stimulus* scanPointStimulus(OBJECT* obj, const Anatomy& anatomy);
   Stimulus* scanTestStimulus(OBJECT* obj);
   Stimulus* scanBoxStimulus(OBJECT* obj, const Anatomy& anatomy);
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
   else if (method == "box")
      return scanBoxStimulus(obj, anatomy);

   assert(false); // reachable only due to bad input
   return 0;
}


namespace
{
   Stimulus* scanPointStimulus(OBJECT* obj, const Anatomy& anatomy)
   {
      PointStimulusParms p;
      objectGet(obj, "cell",   p.cell,   "0");
      objectGet(obj, "freq",   p.freq,   "1000");
      objectGet(obj, "duration",   p.duration,   "1");
      objectGet(obj, "iStim",  p.iStim,  "-52");
      objectGet(obj, "tStart", p.tStart, "0");
      objectGet(obj, "tStop",   p.tStop,   "-1");
      return new PointStimulus(p, anatomy);
   }
}

namespace
{
   Stimulus* scanTestStimulus(OBJECT* obj)
   {
      TestStimulusParms p;
      objectGet(obj, "cell",   p.cell,   "0");
      objectGet(obj, "freq",   p.freq,   "1000");
      objectGet(obj, "iStim",  p.iStim,  "-52");
      objectGet(obj, "rank",   p.rank,   "0");
      objectGet(obj, "tEnd",   p.tEnd,   "2");
      objectGet(obj, "tStart", p.tStart, "1");
      return new TestStimulus(p);
   }
}

namespace
{
   Stimulus* scanBoxStimulus(OBJECT* obj, const Anatomy& anatomy)
   {
      stringstream buf;
      buf << anatomy.nx(); string nxString = buf.str(); buf.str(string());
      buf << anatomy.ny(); string nyString = buf.str(); buf.str(string());
      buf << anatomy.nz(); string nzString = buf.str(); buf.str(string());

      BoxStimulusParms p;
      objectGet(obj, "tStart",   p.tStart,   "1");
      objectGet(obj, "duration", p.duration, "1");
      objectGet(obj, "freq",     p.freq,     "1000");
      objectGet(obj, "iStim",    p.iStim,    "-52");
      objectGet(obj, "xMin",     p.xMin,     "-1");
      objectGet(obj, "yMin",     p.yMin,     "-1");
      objectGet(obj, "zMin",     p.zMin,     "-1");
      objectGet(obj, "xMax",     p.xMax,     nxString);
      objectGet(obj, "yMax",     p.yMax,     nyString);
      objectGet(obj, "zMax",     p.zMax,     nzString);
      
      return new BoxStimulus(p, anatomy);
   }
}
