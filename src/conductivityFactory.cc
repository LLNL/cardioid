#include "conductivityFactory.hh"
#include <cassert>
#include "Conductivity.hh"
#include "object_cc.hh"

#include "FibreConductivity.hh"
#include "UniformConductivity.hh"

using namespace std;

namespace 
{
   Conductivity* scanFibreConductivity(OBJECT* obj);
   Conductivity* scanUniformConductivity(OBJECT* obj);
}


Conductivity* conductivityFactory(const std::string& name)
{
   OBJECT* obj = objectFind(name, "CONDUCTIVITY");
   string method; objectGet(obj, "method", method, "undefined");

   if (method == "undefined")
      assert(1==0);
   else if (method == "fibre")    return scanFibreConductivity(obj);
   else if (method == "fiber")    return scanFibreConductivity(obj);
   else if (method == "uniform")  return scanUniformConductivity(obj);
//   else if (method == "JHU")      return scanJhuConductivity(obj);
   
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
      objectGet(obj, "s11", p.sigma.a11, "0.1");
      objectGet(obj, "s22", p.sigma.a22, "0.1");
      objectGet(obj, "s33", p.sigma.a33, "0.1");
      objectGet(obj, "s12", p.sigma.a12, "0.0");
      objectGet(obj, "s13", p.sigma.a13, "0.0");
      objectGet(obj, "s23", p.sigma.a23, "0.0");
      return new UniformConductivity(p);
   }
}

   
