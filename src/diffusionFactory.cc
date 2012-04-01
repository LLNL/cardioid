#include "diffusionFactory.hh"

#include <cassert>

#include "Diffusion.hh"
#include "object_cc.hh"

#include "Saleheen98Diffusion.hh"
#include "SaleheenDev.hh"
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
                                  const ThreadTeam& threadInfo, int simLoopType);
   Diffusion* saleheen98DiffusionFactory(OBJECT* obj, const Anatomy& anatomy);
   Diffusion* saleheenDevFactory(OBJECT* obj, const Anatomy& anatomy);
   void checkForObsoleteKeywords(OBJECT* obj);
}


Diffusion* diffusionFactory(const string& name, const Anatomy& anatomy,
                            const ThreadTeam& threadInfo, int simLoopType)
{
   OBJECT* obj = objectFind(name, "DIFFUSION");

   checkForObsoleteKeywords(obj);

   string method; objectGet(obj, "method", method, "");

   if (method.empty())
      assert(1==0);
   else if (method == "FGR")
      return fgrDiffusionFactory(obj, anatomy, threadInfo, simLoopType);
   else if (method == "Saleheen98")
      return saleheen98DiffusionFactory(obj, anatomy);
   else if (method == "SaleheenDev")
      return saleheenDevFactory(obj, anatomy);
   else if (method == "null")
      return new NullDiffusion();
   
   assert(false); // reachable only due to bad input
   return 0;
}

namespace
{
   Diffusion* fgrDiffusionFactory(OBJECT* obj, const Anatomy& anatomy,
                                  const ThreadTeam& threadInfo, int simLoopType)
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
         return new FGRDiffusionThreads(p, anatomy, threadInfo);
      if (variant == "simd")
         return new FGRDiffusion(p, anatomy, threadInfo);

      // unreachable.  Should have matched a clause above.
      assert(false);
      return 0; 
   }
}


namespace
{
   Diffusion* saleheen98DiffusionFactory(OBJECT* obj, const Anatomy& anatomy)
   {
      Saleheen98DiffusionParms p;
      objectGet(obj, "diffusionScale", p.diffusionScale_, "1.0", "l^3/capacitance");

      string variant;
      objectGet(obj, "variant", variant, "precompute");

      if (variant == "precompute")
         return new Saleheen98PrecomputeDiffusion(p, anatomy);

      assert(1==0); // reachable only due to bad input.
      return 0;
   }
}

namespace
{
   Diffusion* saleheenDevFactory(OBJECT* obj, const Anatomy& anatomy)
   {
      SaleheenDevParms p;
      objectGet(obj, "diffusionScale", p.diffusionScale_, "1.0", "l^3/capacitance");

      string variant;
      objectGet(obj, "variant", variant, "precompute");

      if (variant == "precompute")
         return new SaleheenDev(p, anatomy);

      assert(1==0); // reachable only due to bad input.
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
