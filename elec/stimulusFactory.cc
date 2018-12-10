#include "stimulusFactory.hh"

#include <cassert>
#include <sstream>
#include "object_cc.hh"
#include "Anatomy.hh"
#include "Stimulus.hh"
#include "PointStimulus.hh"
#include "TestStimulus.hh"
#include "BoxStimulus.hh"
#include "PeriodicPulse.hh"
#include "RandomPulse.hh"

using namespace std;

namespace
{
   Stimulus* scanPointStimulus(OBJECT* obj, const StimulusBaseParms& p, const Anatomy& anatomy, Pulse* pulse);
   Stimulus* scanTestStimulus(OBJECT* obj, const StimulusBaseParms& p, Pulse* pulse);
   Stimulus* scanBoxStimulus(OBJECT* obj, const StimulusBaseParms& p, const Anatomy& anatomy, Pulse* pulse, const std::string& name);
}


Stimulus* stimulusFactory(const std::string& name, const Anatomy& anatomy)
{
   OBJECT* obj = objectFind(name, "STIMULUS");
   string method;
   StimulusBaseParms p;
   objectGet(obj, "method", method, "undefined");
   objectGet(obj, "t0", p.t0, "-1000", "t");
   objectGet(obj, "tf", p.tf, "1e30",  "t");

   string pulse_type;
   double duration;
   double vStim;
   double tStart;
   Pulse* pulse;
   objectGet(obj, "pulse",    pulse_type, "periodic");
   objectGet(obj, "duration", duration, "1",    "t");
   objectGet(obj, "tStart",   tStart,   "0",    "t");
   objectGet(obj, "vStim",    vStim,    "-52",  "voltage/t");
   if (pulse_type == "periodic")
   {
      double period;
      objectGet(obj, "period",   period,   "1000", "t");
      pulse=new PeriodicPulse(period, duration, -vStim, tStart);
   }
   else if (pulse_type == "random")
   {
      double minperiod, maxperiod;
      objectGet(obj, "min_period",   minperiod,   "1000", "t");
      objectGet(obj, "max_period",   maxperiod,   "1000", "t");
      pulse=new RandomPulse(minperiod, maxperiod, duration, -vStim, tStart);
   }

   if (method == "undefined")
      assert(false);
   else if (method == "point")
      return scanPointStimulus(obj, p, anatomy, pulse);
   else if (method == "test")
      return scanTestStimulus(obj, p, pulse);
   else if (method == "box")
      return scanBoxStimulus(obj, p, anatomy, pulse, name);

   assert(false); // reachable only due to bad input
   return 0;
}


namespace
{
   Stimulus* scanPointStimulus(OBJECT* obj, const StimulusBaseParms& bp, const Anatomy& anatomy, Pulse* pulse)
   {
      PointStimulusParms p;
      p.baseParms = bp;
      objectGet(obj, "cell",     p.cell,     "0");
      return new PointStimulus(p, anatomy, pulse);
   }
}

namespace
{
   Stimulus* scanTestStimulus(OBJECT* obj, const StimulusBaseParms& bp, Pulse* pulse)
   {
      TestStimulusParms p;
      p.baseParms = bp;
      objectGet(obj, "cell",     p.cell,     "0");
      objectGet(obj, "rank",     p.rank,     "0");
      return new TestStimulus(p, pulse);
   }
}

namespace
{
   Stimulus* scanBoxStimulus(OBJECT* obj, const StimulusBaseParms& bp, const Anatomy& anatomy, Pulse* pulse, const std::string& name)
   {
      stringstream buf;
      buf << anatomy.nx(); string nxString = buf.str(); buf.str(string());
      buf << anatomy.ny(); string nyString = buf.str(); buf.str(string());
      buf << anatomy.nz(); string nzString = buf.str(); buf.str(string());

      BoxStimulusParms p;
      p.baseParms = bp;
      // TO DO: ELIMINATE GID RANGES
      objectGet(obj, "xMin",     p.xMin,     "-1");
      objectGet(obj, "yMin",     p.yMin,     "-1");
      objectGet(obj, "zMin",     p.zMin,     "-1");
      objectGet(obj, "xMax",     p.xMax,     nxString);
      objectGet(obj, "yMax",     p.yMax,     nyString);
      objectGet(obj, "zMax",     p.zMax,     nzString);

      objectGet(obj, "isXYZ",   p.isXYZ,     "0");
      objectGet(obj, "x",     p.x,     "0");
      objectGet(obj, "y",     p.y,     "0");
      objectGet(obj, "z",     p.z,     "0");
      objectGet(obj, "length",     p.length,     "3.0");

      return new BoxStimulus(p, anatomy, pulse, name);
   }
}
