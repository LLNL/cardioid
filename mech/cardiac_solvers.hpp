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

   // Backtracing line search with (currently disabled) Armijo condition
   virtual double ComputeScalingFactor(const Vector &x, const Vector &b) const;

};

// Custom block preconditioner for the Jacobian of the incompressible nonlinear
// elasticity operator. It has the form
//
// P^-1 = [ K^-1 0 ][ I -B^T ][ I  0           ]
//        [ 0    I ][ 0  I   ][ 0 -\gamma S^-1 ]
//
// where the original Jacobian has the form
//
// J = [ K B^T ]
//     [ B 0   ]
//
// and K^-1 is an approximation of the inverse of the displacement part of the
// Jacobian and S^-1 is an approximation of the inverse of the Schur
// complement S = B K^-1 B^T. The Schur complement is approximiated using
// a mass matrix of the pressure variables.
class JacobianPreconditioner : public Solver
{
protected:
   // Finite element spaces for setting up preconditioner blocks
   Array<ParFiniteElementSpace *> spaces;

   // Offsets for extracting block vector segments
   Array<int> &block_trueOffsets;

   // Jacobian for block access
   BlockOperator *jacobian;

   // Scaling factor for the pressure mass matrix in the block preconditioner
   double gamma;

   // Objects for the block preconditioner application
   Operator *pressure_mass;
   Solver *mass_pcg;
   Solver *mass_prec;
   Solver *stiff_pcg;
   Solver *stiff_prec;

public:
   JacobianPreconditioner(Array<ParFiniteElementSpace *> &fes,
                          Operator &mass, Array<int> &offsets);

   virtual void Mult(const Vector &k, Vector &y) const;
   virtual void SetOperator(const Operator &op);

   virtual ~JacobianPreconditioner();
};
   
}

#endif
