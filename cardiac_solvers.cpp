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

   add(x, -scale, c, x_out);
   oper->Mult(x_out, test);

   double norm = Norm(test);

   int myid;
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   
   // Decreases the scaling of the update until the new mesh is valid.
   while (norm > 1.2 * norm0 || isnan(norm) != 0) {
      scale *= factor;
      add(x, -scale, c, x_out);
      oper->Mult(x_out, test);
      norm = Norm(test);
   }
   return scale;
}

}
