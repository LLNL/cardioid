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

   /// Armijo constant
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

JacobianPreconditioner::JacobianPreconditioner(Array<ParFiniteElementSpace *>
                                               &fes,
                                               Operator &mass,
                                               Array<int> &offsets)
   : Solver(offsets[2]), block_trueOffsets(offsets), pressure_mass(&mass)
{
   fes.Copy(spaces);

   gamma = 0.00001;

   // The mass matrix and preconditioner do not change every Newton cycle, so
   // we only need to define them once
   HypreBoomerAMG *mass_prec_amg = new HypreBoomerAMG();
   mass_prec_amg->SetPrintLevel(0);

   mass_prec = mass_prec_amg;

   CGSolver *mass_pcg_iter = new CGSolver(spaces[0]->GetComm());
   mass_pcg_iter->SetRelTol(1e-12);
   mass_pcg_iter->SetAbsTol(1e-12);
   mass_pcg_iter->SetMaxIter(200);
   mass_pcg_iter->SetPrintLevel(0);
   mass_pcg_iter->SetPreconditioner(*mass_prec);
   mass_pcg_iter->SetOperator(*pressure_mass);
   mass_pcg_iter->iterative_mode = false;

   mass_pcg = mass_pcg_iter;

   // The stiffness matrix does change every Newton cycle, so we will define it
   // during SetOperator
   stiff_pcg = NULL;
   stiff_prec = NULL;
}

void JacobianPreconditioner::Mult(const Vector &k, Vector &y) const
{
   // Extract the blocks from the input and output vectors
   Vector disp_in(k.GetData() + block_trueOffsets[0],
                  block_trueOffsets[1]-block_trueOffsets[0]);
   Vector pres_in(k.GetData() + block_trueOffsets[1],
                  block_trueOffsets[2]-block_trueOffsets[1]);

   Vector disp_out(y.GetData() + block_trueOffsets[0],
                   block_trueOffsets[1]-block_trueOffsets[0]);
   Vector pres_out(y.GetData() + block_trueOffsets[1],
                   block_trueOffsets[2]-block_trueOffsets[1]);

   Vector temp(block_trueOffsets[1]-block_trueOffsets[0]);
   Vector temp2(block_trueOffsets[1]-block_trueOffsets[0]);

   // Perform the block elimination for the preconditioner
   mass_pcg->Mult(pres_in, pres_out);
   pres_out *= -gamma;

   jacobian->GetBlock(0,1).Mult(pres_out, temp);
   subtract(disp_in, temp, temp2);

   stiff_pcg->Mult(temp2, disp_out);
}

void JacobianPreconditioner::SetOperator(const Operator &op)
{
   jacobian = (BlockOperator *) &op;

   // Initialize the stiffness preconditioner and solver
   if (stiff_prec == NULL)
   {
      HypreBoomerAMG *stiff_prec_amg = new HypreBoomerAMG();
      stiff_prec_amg->SetPrintLevel(0);
      stiff_prec_amg->SetElasticityOptions(spaces[0]);

      stiff_prec = stiff_prec_amg;

      GMRESSolver *stiff_pcg_iter = new GMRESSolver(spaces[0]->GetComm());
      stiff_pcg_iter->SetRelTol(1e-4);
      stiff_pcg_iter->SetAbsTol(1e-4);
      stiff_pcg_iter->SetMaxIter(2000);
      stiff_pcg_iter->SetPrintLevel(0);
      stiff_pcg_iter->SetPreconditioner(*stiff_prec);
      stiff_pcg_iter->iterative_mode = false;

      stiff_pcg = stiff_pcg_iter;
   }

   // At each Newton cycle, compute the new stiffness AMG preconditioner by
   // updating the iterative solver which, in turn, updates its preconditioner
   stiff_pcg->SetOperator(jacobian->GetBlock(0,0));
}

JacobianPreconditioner::~JacobianPreconditioner()
{
   delete mass_pcg;
   delete mass_prec;
   delete stiff_prec;
   delete stiff_pcg;
}
   
}
