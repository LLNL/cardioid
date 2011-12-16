#include "diffusionFactory.hh"

#include <cassert>

#include "Diffusion.hh"
#include "object_cc.hh"


#include "Saleheen98Diffusion.hh"

class Anatomy;



using namespace std;


namespace
{
   Diffusion* saleheen98DiffusionFactory(OBJECT* obj, const Anatomy& anatomy);
}


Diffusion* diffusionFactory(const string& name, const Anatomy& anatomy)
{
   OBJECT* obj = objectFind(name, "DIFFUSION");
   string method; objectGet(obj, "method", method, "undefined");

   if (method == "undefined")
      assert(1==0);
   else if (method == "Saleheen98")
      return saleheen98DiffusionFactory(obj, anatomy);
   
   assert(false); // reachable only due to bad input
   return 0;
}



namespace
{
   Diffusion* saleheen98DiffusionFactory(OBJECT* obj, const Anatomy& anatomy)
   {
      Saleheen98DiffusionParms p;
      objectGet(obj, "conductivity", p.conductivityName_, "conductivity");
      objectGet(obj, "diffusionScale", p.diffusionScale_, "1.0");

      string variant;
      objectGet(obj, "variant", variant, "precompute");

      if (variant == "precompute")
         return new Saleheen98PrecomputeDiffusion(p, anatomy);

      assert(1==0); // reachable only due to bad input.
      return 0;
   }
}

