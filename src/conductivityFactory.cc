#include "conductivityFactory.hh"
#include <cassert>
#include "Conductivity.hh"
#include "object_cc.hh"

#include "Anatomy.hh"
#include "FibreConductivity.hh"
#include "UniformConductivity.hh"
#include "JHUConductivity.hh"

using namespace std;

namespace 
{
   Conductivity* scanFibreConductivity(OBJECT* obj);
   Conductivity* scanUniformConductivity(OBJECT* obj);
   Conductivity* scanJHUConductivity(OBJECT* obj, const Anatomy& anatomy);
}


Conductivity* conductivityFactory(const string& name, const Anatomy& anatomy)
{
   OBJECT* obj = objectFind(name, "CONDUCTIVITY");
   string method; objectGet(obj, "method", method, "undefined");

   if (method == "undefined")
      assert(1==0);
   else if (method == "fibre")    return scanFibreConductivity(obj);
   else if (method == "fiber")    return scanFibreConductivity(obj);
   else if (method == "uniform")  return scanUniformConductivity(obj);
   else if (method == "JHU")      return scanJHUConductivity(obj, anatomy);
   
   assert(1==0); // reachable only due to bad input

}


namespace
{
   Conductivity* scanFibreConductivity(OBJECT* obj)
   {
      FibreConductivityParms p;
      objectGet(obj, "sigmaTi", p.sigmaTi, "0.0315");
      objectGet(obj, "sigmaLi", p.sigmaLi, "0.3");
      return new FibreConductivity(p);
   }
}

namespace
{
   Conductivity* scanUniformConductivity(OBJECT* obj)
   {
      UniformConductivityParms p;
      objectGet(obj, "sigma11", p.sigma.a11, "0.1");
      objectGet(obj, "sigma22", p.sigma.a22, "0.1");
      objectGet(obj, "sigma33", p.sigma.a33, "0.1");
      objectGet(obj, "sigma12", p.sigma.a12, "0.0");
      objectGet(obj, "sigma13", p.sigma.a13, "0.0");
      objectGet(obj, "sigma23", p.sigma.a23, "0.0");
      return new UniformConductivity(p);
   }
}

namespace
{
   /** See JHUConductivity.hh for an explanation of the paramters and
    *  their default values. */
   Conductivity* scanJHUConductivity(OBJECT* obj, const Anatomy& anatomy)
   {
      JHUConductivityParms p;
      p.nx = anatomy.nx();
      p.ny = anatomy.ny();
      p.nz = anatomy.nz();
      p.transmuralAxis = anatomy.ny() - 2;

      objectGet(obj, "sigmaTi", p.sigmaTi, "0.0315");
      objectGet(obj, "sigmaLi", p.sigmaLi, "0.3");

      objectGet(obj, "sheetAngle", p.sheetAngle, "45");
      objectGet(obj, "rotatationMatrix", p.rotatationMatrix, "1");
      objectGet(obj, "homogeneousFiber", p.homogeneousFiber, "0");

      return new JHUConductivity(p);
   }
}

