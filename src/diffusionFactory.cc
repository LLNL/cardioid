#include "diffusionFactory.hh"

#include <cassert>

#include "Diffusion.hh"
#include "object_cc.hh"


#include "Salheen98Diffusion.hh"

class Anatomy;



using namespace std;


namespace
{
   Diffusion* salheen98DiffusionFactory(OBJECT* obj, const Anatomy& anatomy);
}


Diffusion* diffusionFactory(const string& name, const Anatomy& anatomy)
{
   OBJECT* obj = objectFind(name, "DIFFUSION");
   string method; objectGet(obj, "method", method, "undefined");

   if (method == "undefined")
      assert(1==0);
   else if (method == "Salheen98")
      return salheen98DiffusionFactory(obj, anatomy);
   
   assert(1==0); // reachable only due to bad input
}



namespace
{
   Diffusion* salheen98DiffusionFactory(OBJECT* obj, const Anatomy& anatomy)
   {
      Salheen98DiffusionParms p;
      objectGet(obj, "conductivity", p.conductivityName_, "conductivity");

      string variant;
      objectGet(obj, "variant", variant, "precompute");

      if (variant == "precompute")
	 return new Salheen98PrecomputeDiffusion(p, anatomy);

      assert(1==0); // reachable only due to bad input.
   }
}

