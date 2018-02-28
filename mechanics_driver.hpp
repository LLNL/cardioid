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

   ActiveTensionFunction *tension_func;
   QuadratureSpace *Q_space;
   
   ParBlockNonlinearForm *Hform;

   /// Newton solver for the hyperelastic operator
   CardiacNewtonSolver newton_solver;
   /// Solver for the Jacobian solve in the Newton method
   Solver *J_solver;
   /// Preconditioner for the Jacobian
   Solver *J_prec;
   /// Specific cardiac hyperelastic model
   CardiacModel *model;

   bool slu_solver;

   VectorFunctionCoefficient *fib;
   FunctionCoefficient *pres;
   VectorFunctionCoefficient *vol;
   QuadratureFunctionCoefficient *qat;

   double dt;
   
public:
   CardiacOperator(Array<ParFiniteElementSpace *> &fes, Array<Array<int>*> &ess_bdr, Array<int> &pres_bdr, bool slu, Array<int> &block_trueOffsets, double rel_tol, double abs_tol, int iter, double timestep);

   /// Required to use the native newton solver
   virtual Operator &GetGradient(const Vector &xp) const;
   virtual void Mult(const Vector &k, Vector &y) const;

   /// Driver for the newton solver
   void Solve(Vector &xp) const;
   
   virtual ~CardiacOperator();
};

void visualize(ostream &out, ParMesh *mesh, ParGridFunction *deformed_nodes,
               ParGridFunction *field, const char *field_name = NULL,
               bool init_vis = false);
#endif
