#ifndef MECHANICS_DRIVER
#define MECHANICS_DRIVER

#include "mfem.hpp"
#include "cardiac_physics.hpp"
#include "cardiac_integrators.hpp"
#include "cardiac_solvers.hpp"
#include "cardiac_coefficients.hpp"

using namespace std;
using namespace mfem;


/** After spatial discretization, the cardiac model can be written as:
 *     0=H(x)
 *  where x is the block vector representing the deformation and pressure
 *  and H(x) is the nonlinear cardiac operator. */
class CardiacOperator : public TimeDependentOperator
{
protected:
   Array<ParFiniteElementSpace *> spaces;

   /// Active tension function (using melodee code)
   ActiveTensionFunction *tension_func;

   /// Quadrature space information
   QuadratureSpace *Q_space;

   /// Nonlinear form operator
   ParBlockNonlinearForm *Hform;

   /// Pressure mass for the preconditioner
   Operator *pressure_mass;
   
   /// Newton solver for the hyperelastic operator
   CardiacNewtonSolver newton_solver;

   /// Solver for the Jacobian solve in the Newton method
   Solver *J_solver;

   /// Preconditioner for the Jacobian (not currently used)
   Solver *J_prec;

   /// Specific cardiac hyperelastic model
   CardiacModel *model;

   /// Function coefficients
   VectorFunctionCoefficient *fib;
   FunctionCoefficient *pres;
   VectorFunctionCoefficient *vol;
   QuadratureFunctionCoefficient *qat;
   
   /// Timestep
   double dt;

   /// Direct solver flag
   bool slu;
   
public:
   CardiacOperator(Array<ParFiniteElementSpace *> &fes, Array<Array<int>*> &ess_bdr, Array<int> &pres_bdr, Array<int> &block_trueOffsets, double rel_tol, double abs_tol, int iter, double timestep, bool superlu);

   /// Required to use the native newton solver
   /// Returns the Jacobian matrix (gradient of the residual vector)
   virtual Operator &GetGradient(const Vector &xp) const;
   /// Computes the residual vector
   virtual void Mult(const Vector &k, Vector &y) const;

   /// Driver for the newton solver
   void Solve(Vector &xp) const;
   
   virtual ~CardiacOperator();
};

#endif
