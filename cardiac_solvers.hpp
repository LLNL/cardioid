#ifndef CARDIAC_SOLVE
#define CARDIAC_SOLVE

#include "mfem.hpp"

namespace mfem
{

/// Enhanced Cardiac Newton solver
class CardiacNewtonSolver : public NewtonSolver
{
private:
   // line search scaling factor
   const double factor;
   
public:
   CardiacNewtonSolver(MPI_Comm _comm, double _fac = 0.5)
      : NewtonSolver(_comm), factor(_fac) { }
   
   virtual double ComputeScalingFactor(const Vector &x, const Vector &b) const;

};

}

#endif
