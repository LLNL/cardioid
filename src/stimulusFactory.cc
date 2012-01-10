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
   Stimulus* scanPointStimulus(OBJECT* obj, const StimulusBaseParms& p, const Anatomy& anatomy);
   Stimulus* scanTestStimulus(OBJECT* obj, const StimulusBaseParms& p);
   Stimulus* scanBoxStimulus(OBJECT* obj, const StimulusBaseParms& p, const Anatomy& anatomy);
}


Stimulus* stimulusFactory(const std::string& name, const Anatomy& anatomy)
{
   OBJECT* obj = objectFind(name, "STIMULUS");
   string method;
   StimulusBaseParms p;
   objectGet(obj, "method", method, "undefined");
   objectGet(obj, "t0", p.t0, "-1000");
   objectGet(obj, "tf", p.tf, "1e30");
   
   if (method == "undefined")
      assert(false);
   else if (method == "point")
      return scanPointStimulus(obj, p, anatomy);
   else if (method == "test")
      return scanTestStimulus(obj, p);
   else if (method == "box")
      return scanBoxStimulus(obj, p, anatomy);

   assert(false); // reachable only due to bad input
   return 0;
}


namespace
{
   Stimulus* scanPointStimulus(OBJECT* obj, const StimulusBaseParms& bp, const Anatomy& anatomy)
   {
      PointStimulusParms p;
      p.baseParms = bp;
      objectGet(obj, "cell",     p.cell,     "0");
      objectGet(obj, "duration", p.duration, "1");
      objectGet(obj, "period",   p.period,   "1000");
      objectGet(obj, "tStart",   p.tStart,   "0");
      objectGet(obj, "vStim",    p.vStim,    "-52");
      return new PointStimulus(p, anatomy);
   }
}

namespace
{
   Stimulus* scanTestStimulus(OBJECT* obj, const StimulusBaseParms& bp)
   {
      TestStimulusParms p;
      p.baseParms = bp;
      objectGet(obj, "cell",     p.cell,     "0");
      objectGet(obj, "rank",     p.rank,     "0");
      objectGet(obj, "duration", p.duration, "1");
      objectGet(obj, "period",   p.period,   "1000");
      objectGet(obj, "tStart",   p.tStart,   "1");
      objectGet(obj, "vStim",    p.vStim,    "-52");
      return new TestStimulus(p);
   }
}

namespace
{
   Stimulus* scanBoxStimulus(OBJECT* obj, const StimulusBaseParms& bp, const Anatomy& anatomy)
   {
      stringstream buf;
      buf << anatomy.nx(); string nxString = buf.str(); buf.str(string());
      buf << anatomy.ny(); string nyString = buf.str(); buf.str(string());
      buf << anatomy.nz(); string nzString = buf.str(); buf.str(string());

      BoxStimulusParms p;
      p.baseParms = bp;
      objectGet(obj, "duration", p.duration, "1");
      objectGet(obj, "period",   p.period,   "1000");
      objectGet(obj, "tStart",   p.tStart,   "1");
      objectGet(obj, "vStim",    p.vStim,    "-52");
      objectGet(obj, "xMin",     p.xMin,     "-1");
      objectGet(obj, "yMin",     p.yMin,     "-1");
      objectGet(obj, "zMin",     p.zMin,     "-1");
      objectGet(obj, "xMax",     p.xMax,     nxString);
      objectGet(obj, "yMax",     p.yMax,     nyString);
      objectGet(obj, "zMax",     p.zMax,     nzString);
      
      return new BoxStimulus(p, anatomy);
   }
}
