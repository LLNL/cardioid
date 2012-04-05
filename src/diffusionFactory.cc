#include "diffusionFactory.hh"

#include <cassert>

#include "Diffusion.hh"
#include "object_cc.hh"

#include "FGRDiffusion.hh"
#include "FGRDiffusionOMP.hh"
#include "FGRDiffusionThreads.hh"
#include "NullDiffusion.hh"

class Anatomy;
class ThreadTeam;


using namespace std;


namespace
{
   Diffusion* fgrDiffusionFactory(OBJECT* obj, const Anatomy& anatomy,
                                  const ThreadTeam& threadInfo,
                                  const ThreadTeam& reactionThreadInfo,
                                  int simLoopType);
   void checkForObsoleteKeywords(OBJECT* obj);
}


Diffusion* diffusionFactory(const string& name, const Anatomy& anatomy,
                            const ThreadTeam& threadInfo,
                            const ThreadTeam& reactionThreadInfo,
                            int simLoopType)
{
   OBJECT* obj = objectFind(name, "DIFFUSION");

   checkForObsoleteKeywords(obj);

   string method; objectGet(obj, "method", method, "FGR");

   if (method.empty())
      assert(1==0);
   else if (method == "FGR")
      return fgrDiffusionFactory(obj, anatomy, threadInfo, reactionThreadInfo, simLoopType);
   else if (method == "null")
      return new NullDiffusion();
   
   assert(false); // reachable only due to bad input
   return 0;
}

namespace
{
   Diffusion* fgrDiffusionFactory(OBJECT* obj, const Anatomy& anatomy,
                                  const ThreadTeam& threadInfo,
                                  const ThreadTeam& reactionThreadInfo,
                                  int simLoopType)
   {
      FGRUtils::FGRDiffusionParms p;
      objectGet(obj, "diffusionScale", p.diffusionScale_, "1.0", "l^3/capacitance");
      string defaultVariant = "omp";
      if (simLoopType == 1) // parallelDiffusionReaction
         defaultVariant = "simd";
      string variant;
      objectGet(obj, "variant", variant, defaultVariant);
      if (variant == "omp")
         return new FGRDiffusionOMP(p, anatomy);
      if (variant == "threads")
         return new FGRDiffusionThreads(p, anatomy, threadInfo, reactionThreadInfo);
      if (variant == "simd")
         return new FGRDiffusion(p, anatomy, threadInfo, reactionThreadInfo);

      // unreachable.  Should have matched a clause above.
      std::cerr<<"ERROR --- invalid 'variant' parameter: "<<variant<<std::endl;
      assert(false);
      return 0; 
   }
}


namespace
{
   void checkForObsoleteKeywords(OBJECT* obj)
   {
      assert (object_testforkeyword(obj, "conductivity") == 0);
   }
}
