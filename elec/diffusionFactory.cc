#include "diffusionFactory.hh"

#include <cassert>

#include "Diffusion.hh"
#include "object_cc.hh"

#include "FGRUtils.hh"
#include "FGRDiffusion.hh"
#include "FGRDiffusionOMP.hh"
#include "FGRDiffusionThreads.hh"
#include "FGRDiffusionStrip.hh"
#include "FGRDiffusionOverlap.hh"
#include "NullDiffusion.hh"
//#include "OpenmpGpuRedblackDiffusion.hh"
//#include "OpenmpGpuFlatDiffusion.hh"
#include "CUDADiffusion.hh"
#include "Simulate.hh"

class Anatomy;
class ThreadTeam;


using namespace std;


namespace
{
   Diffusion* fgrDiffusionFactory(OBJECT* obj, const Anatomy& anatomy,
                                  const ThreadTeam& threadInfo,
                                  const ThreadTeam& reactionThreadInfo,
                                  int simLoopType, string loadLevelVariant);
   void checkForObsoleteKeywords(OBJECT* obj);
}


Diffusion* diffusionFactory(const string& name, const Anatomy& anatomy,
                            const ThreadTeam& threadInfo,
                            const ThreadTeam& reactionThreadInfo,
                            int simLoopType, string &variantHint)
{
   OBJECT* obj = objectFind(name, "DIFFUSION");

   checkForObsoleteKeywords(obj);

   string method; objectGet(obj, "method", method, "FGR");
   double diffusionScale;
   objectGet(obj, "diffusionScale", diffusionScale, "1.0", "l^3/capacitance");
      
   if (method.empty())
      assert(1==0);
   else if (method == "FGR")
      return fgrDiffusionFactory(obj, anatomy, threadInfo, reactionThreadInfo, simLoopType, variantHint);
   //else if (method == "gpu" || method == "OpenmpGpuRedblack")
   //   return new OpenmpGpuRedblackDiffusion(anatomy, simLoopType);
   //else if (method == "OpenmpGpuFlat")
   //   return new OpenmpGpuFlatDiffusion(anatomy, simLoopType);
#ifdef USE_CUDA
   else if (method == "cuda")
      return new CUDADiffusion(anatomy, simLoopType, diffusionScale); 
#endif
   else if (method == "null")
      return new NullDiffusion(anatomy, simLoopType);
   
   assert(false); // reachable only due to bad input
   return 0;
}

namespace
{
   Diffusion* fgrDiffusionFactory(OBJECT* obj, const Anatomy& anatomy,
                                  const ThreadTeam& threadInfo,
                                  const ThreadTeam& reactionThreadInfo,
                                  int simLoopType, string variantHint)
   {
      FGRUtils::FGRDiffusionParms p;
      objectGet(obj, "diffusionScale", p.diffusionScale_, "1.0", "l^3/capacitance");
      objectGet(obj, "printBBox",      p.printBBox_, "0");
      string defaultVariant = "omp";
      if (! variantHint.empty())
         defaultVariant = variantHint;
      string variant;
      objectGet(obj, "variant", variant, defaultVariant);
      if (0) {}
      else if (variant == "omp" || simLoopType == Simulate::omp)
         return new FGRDiffusionOMP(p, anatomy);
      else if (variant == "threads")
         return new FGRDiffusionThreads(p, anatomy, threadInfo, reactionThreadInfo);
      else if (variant == "simd")
         return new FGRDiffusion(p, anatomy, threadInfo, reactionThreadInfo);
      else if (variant == "strip" )
         return new FGRDiffusionStrip(p, anatomy, threadInfo, reactionThreadInfo);
      else if (variant == "overlap" )
         return new FGRDiffusionOverlap(p, anatomy, threadInfo, reactionThreadInfo);


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
