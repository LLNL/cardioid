#include "mfem.hpp"
#include "cardiac_solvers.hpp"

namespace mfem
{

double CardiacNewtonSolver::ComputeScalingFactor(const Vector &x,
                                                 const Vector &b) const
{
   
   Vector x_out(x.Size());
   Vector test(x.Size());
   double scale = 1.0;
   double norm0 = Norm(r);

   double m = Dot(c,c);

   double arm_c = 0.0;
   
   add(x, -scale, c, x_out);
   oper->Mult(x_out, test);

   double norm = Norm(test);

   
   while (isnan(norm) != 0 || (norm > 1.2 * (norm0 - scale * arm_c * m) && scale > 0.001)) {

      scale *= factor;
      add(x, -scale, c, x_out);
      oper->Mult(x_out, test);
      norm = Norm(test);
      
   }
   

   return scale;
}

}
