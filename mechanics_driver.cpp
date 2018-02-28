#include "mfem.hpp"
#include "cardiac_physics.hpp"
#include "cardiac_integrators.hpp"
#include "mechanics_driver.hpp"
#include "ConstantTension.hpp"
#include "Lumens2009.hpp"
#include <memory>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <cerrno>

int main(int argc, char *argv[])
{
   // Initialize MPI.
   int num_procs, myid;
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);

   // Parse command-line options.
   const char *mesh_file = "./beam-hex.mesh";
   int ser_ref_levels = 0;
   int par_ref_levels = 0;
   int order = 2;
   bool slu_solver = true;
   bool visualization = true;
   double newton_rel_tol = 1.0e-2;
   double newton_abs_tol = 1.0e-12;
   int newton_iter = 500;
   double tf = 1.0;
   double dt = 1.0;

   OptionsParser args(argc, argv);
   args.AddOption(&run_mode, "-rm", "--run-mode",
                  "Run mode. 1 is cantilever beam, 2 is inflated ventricle, and 3 is active tension beam.");   
   args.AddOption(&ser_ref_levels, "-rs", "--refine-serial",
                  "Number of times to refine the mesh uniformly in serial.");
   args.AddOption(&par_ref_levels, "-rp", "--refine-parallel",
                  "Number of times to refine the mesh uniformly in parallel.");
   args.AddOption(&order, "-o", "--order",
                  "Order (degree) of the finite elements.");
   args.AddOption(&slu_solver, "-slu", "--superlu", "-no-slu",
                  "--no-superlu", "Use the SuperLU Solver.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&newton_rel_tol, "-rel", "--relative-tolerance",
                  "Relative tolerance for the Newton solve.");
   args.AddOption(&newton_abs_tol, "-abs", "--absolute-tolerance",
                  "Absolute tolerance for the Newton solve.");
   args.AddOption(&newton_iter, "-it", "--newton-iterations",
                  "Maximum iterations for the Newton solve.");
   args.AddOption(&tf, "-tf", "--t-final",
                  "Final time.");
   args.AddOption(&dt, "-dt", "--time-step",
                  "Length of time step.");
   
   args.Parse();
   if (!args.Good())
   {
      if (myid == 0)
      {
         args.PrintUsage(cout);
      }
      MPI_Finalize();
      return 1;
   }
   if (myid == 0)
   {
      args.PrintOptions(cout);
   }
   
   // Open the mesh
   Mesh *mesh = NULL;
   if (run_mode == 1 || run_mode == 3) {
      const char *mesh_file = "./beam-hex.mesh";
      mesh = new Mesh(mesh_file, 1, 1);
   }
   if (run_mode == 2) {
      const char *mesh_file = "./hollow-ball.vtk";
      mesh = new Mesh(mesh_file, 1, 1);
      setSurfaces(mesh);
   }

   ParMesh *pmesh = NULL;

   for (int lev = 0; lev < ser_ref_levels; lev++)
      {
         mesh->UniformRefinement();
      }
   pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   for (int lev = 0; lev < par_ref_levels; lev++)
      {
         pmesh->UniformRefinement();
      }

   delete mesh;
   int dim = pmesh->Dimension();

   
   // Define the finite element spaces for displacement and pressure (Stokes elements)
   H1_FECollection quad_coll(order, dim);
   H1_FECollection lin_coll(order-1, dim);
   ParFiniteElementSpace R_space(pmesh, &quad_coll, dim);
   ParFiniteElementSpace W_space(pmesh, &lin_coll);

   Array<ParFiniteElementSpace *> spaces(2);
   spaces[0] = &R_space;
   spaces[1] = &W_space;

   HYPRE_Int glob_R_size = R_space.GlobalTrueVSize();
   HYPRE_Int glob_W_size = W_space.GlobalTrueVSize();

   // Define the Dirichlet conditions (set to boundary attribute 1)
   Array<Array<int> *> ess_bdr(2);

   // Boundary condition markers
   // 1 - Dirichlet/Essential
   // 2 - Free surface
   // 3 - Traction
   // 4 - Pressure
   Array<int> ess_bdr_u(R_space.GetMesh()->bdr_attributes.Max());
   Array<int> ess_bdr_p(W_space.GetMesh()->bdr_attributes.Max());
   Array<int> pres_bdr(R_space.GetMesh()->bdr_attributes.Max());
     
   ess_bdr_p = 0;
   ess_bdr_u = 0;
   ess_bdr_u[0] = 1;

   pres_bdr = 0;
   pres_bdr[3] = 1;

   ess_bdr[0] = &ess_bdr_u;
   ess_bdr[1] = &ess_bdr_p;


   // Print the mesh statistics
   if (myid == 0)
   {
      std::cout << "***********************************************************\n";
      std::cout << "dim(u) = " << glob_R_size << "\n";
      std::cout << "dim(p) = " << glob_W_size << "\n";
      std::cout << "dim(u+p) = " << glob_R_size + glob_W_size << "\n";
      std::cout << "***********************************************************\n";
   }

   // Define the block structure of the solution vector (u then p)
   Array<int> block_trueOffsets(3); 
   block_trueOffsets[0] = 0;
   block_trueOffsets[1] = R_space.TrueVSize();
   block_trueOffsets[2] = W_space.TrueVSize();
   block_trueOffsets.PartialSum();

   BlockVector xp(block_trueOffsets);

   // Define grid functions for the current configuration, reference configuration,
   // final deformation, and pressure
   ParGridFunction x_gf(&R_space);
   ParGridFunction x_ref(&R_space);
   ParGridFunction x_def(&R_space);
   ParGridFunction p_gf(&W_space);

   // Project the initial and reference configuration functions onto the appropriate grid functions
   VectorFunctionCoefficient deform(dim, InitialDeformation);
   VectorFunctionCoefficient refconfig(dim, ReferenceConfiguration);
  
   x_gf.ProjectCoefficient(deform);
   x_ref.ProjectCoefficient(refconfig);
   p_gf = 0.0;
   
   // Set up the block solution vectors
   x_gf.GetTrueDofs(xp.GetBlock(0));
   p_gf.GetTrueDofs(xp.GetBlock(1));

   // Initialize the cardiac mechanics operator
   CardiacOperator oper(spaces, ess_bdr, pres_bdr, slu_solver, block_trueOffsets,
                        newton_rel_tol, newton_abs_tol, newton_iter, dt);

   // Loop over the timesteps
   for (double t = 0.0; t<tf; t += dt) {

      // Solve the Newton system       
      oper.SetTime(t + dt);
      oper.Solve(xp);
      if (myid == 0) {
         std::cout << "***********************************************************\n";         
         std::cout << "Solve at time " << t + dt << " complete\n";
         std::cout << "***********************************************************\n";         
      }
   }
      
   // Distribute the ghost dofs
   x_gf.Distribute(xp.GetBlock(0));
   p_gf.Distribute(xp.GetBlock(1));

   // Set the final deformation
   subtract(x_gf, x_ref, x_def);

   // Visualize the results if requested
   socketstream vis_u, vis_p;
   if (visualization) {
      char vishost[] = "localhost";
      int  visport   = 19916;
      vis_u.open(vishost, visport);
      vis_u.precision(8);
      visualize(vis_u, pmesh, &x_gf, &x_def, "Deformation", true);
      // Make sure all ranks have sent their 'u' solution before initiating
      // another set of GLVis connections (one from each rank):
      MPI_Barrier(pmesh->GetComm());
      vis_p.open(vishost, visport);
      vis_p.precision(8);
      visualize(vis_p, pmesh, &x_gf, &p_gf, "Pressure", true);
   }

   // Save the displaced mesh, the final deformation, and the pressure
   {
      GridFunction *nodes = &x_gf;
      int owns_nodes = 0;
      pmesh->SwapNodes(nodes, owns_nodes);

      if (myid == 0)
      {
         int err = mkdir("output", 0777);
         err = (err && (errno != EEXIST)) ? 1 : 0;

         if (err == 1) {
            std::cout << "Error creating output directory!\n";
         }         
      }
      
      ostringstream mesh_name, pressure_name, deformation_name;
      mesh_name << "output/mesh." << setfill('0') << setw(6) << myid;
      pressure_name << "output/pressure." << setfill('0') << setw(6) << myid;
      deformation_name << "output/deformation." << setfill('0') << setw(6) << myid;

      ofstream mesh_ofs(mesh_name.str().c_str());
      mesh_ofs.precision(8);
      pmesh->Print(mesh_ofs);
    
      ofstream pressure_ofs(pressure_name.str().c_str());
      pressure_ofs.precision(8);
      p_gf.Save(pressure_ofs);

      ofstream deformation_ofs(deformation_name.str().c_str());
      deformation_ofs.precision(8);
      x_def.Save(deformation_ofs);
   }


   // Free the used memory.
   delete pmesh;
   
   MPI_Finalize();

   return 0;
}

CardiacOperator::CardiacOperator(Array<ParFiniteElementSpace *>&fes,
                                 Array<Array<int> *>&ess_bdr,
                                 Array<int> &pres_bdr,
                                 bool slu,
                                 Array<int> &block_trueOffsets,
                                 double rel_tol,
                                 double abs_tol,
                                 int iter,
                                 double timestep)
   : TimeDependentOperator(fes[0]->TrueVSize() + fes[1]->TrueVSize(), 0.0), 
     newton_solver(fes[0]->GetComm(), 0.8), slu_solver(slu), dt(timestep)
{
   Array<Vector *> rhs(2);
   rhs = NULL;
   tension_func = NULL;
   qat = NULL;
   
   fes.Copy(spaces);

   // Define the quadrature space for the active tension coefficient
   Q_space = new QuadratureSpace(fes[0]->GetMesh(), 2*(fes[0]->GetOrder(0))+3);

   // Define and initialize the quadrature function for the active tension coefficient
   if (run_mode == 3) {
      tension_func = new ActiveTensionFunction(Q_space, fes[0]);
      tension_func->Initialize();
      (*tension_func) = 1.0;
   }
   
   // Define the forcing function coefficients
   fib = new VectorFunctionCoefficient(3, FiberFunction);
   pres = new FunctionCoefficient(PressureFunction);
   vol = new VectorFunctionCoefficient(3, VolumeFunction);

   if (run_mode == 3) {
      qat = new QuadratureFunctionCoefficient(tension_func);
   }
      
   // Initialize the Cardiac model (transversely isotropic)

   if (run_mode == 1) {
      model = new CardiacModel (2.0, 8.0, 2.0, 2.0, 4.0, 4.0, 2.0);
   }
   else {
      model = new CardiacModel (10, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
   }

   // Define the mixed nonlinear form 
   Hform = new ParBlockNonlinearForm(spaces);

   // Add the passive stress integrator
   Hform->AddDomainIntegrator(new CardiacNLFIntegrator(model, *fib));

   // Add the active tension integrator (only mode 3)
   if (run_mode == 3) {
      Hform->AddDomainIntegrator(new ActiveTensionNLFIntegrator(*qat, *fib));
   }
      
   // Add the pressure and traction boundary integrators
   if (run_mode == 1 || run_mode == 2) {
      Hform->AddBdrFaceIntegrator(new PressureBoundaryNLFIntegrator(*pres, *vol), pres_bdr);
   }
   // Set the essential boundary conditions
   Hform->SetEssentialBC(ess_bdr, rhs);

   // Set the solver parameters (only dsuperlu works currently)
   if (slu_solver) {
      SuperLUSolver *superlu = NULL;
      superlu = new SuperLUSolver(MPI_COMM_WORLD);
      superlu->SetPrintStatistics(false);
      superlu->SetSymmetricPattern(false);
      superlu->SetColumnPermutation(superlu::PARMETIS);
      
      J_solver = superlu;
      J_prec = NULL;
   }
   else {
      HypreIdentity *invU, *invP;
      invU = new HypreIdentity;
      invP = new HypreIdentity;

      MINRESSolver *J_minres = new MINRESSolver(MPI_COMM_WORLD);
      J_minres->SetRelTol(rel_tol);
      J_minres->SetAbsTol(0.0);
      J_minres->SetMaxIter(300);
      J_minres->SetPrintLevel(-1);
      J_minres->SetPreconditioner(*J_prec);
      J_solver = J_minres;  
   }

   // Set the newton solve parameters
   newton_solver.iterative_mode = true;
   newton_solver.SetSolver(*J_solver);
   newton_solver.SetOperator(*this);
   newton_solver.SetPrintLevel(1); 
   newton_solver.SetRelTol(rel_tol);
   newton_solver.SetAbsTol(abs_tol);
   newton_solver.SetMaxIter(iter);
}

// Solve the Newton system
void CardiacOperator::Solve(Vector &xp) const
{
   Vector zero;
   newton_solver.Mult(zero, xp);
   if (run_mode == 3) {
      tension_func->CommitStep(dt);
   }
   MFEM_VERIFY(newton_solver.GetConverged(), "Newton Solver did not converge.");
}

// compute: y = H(x,p)
void CardiacOperator::Mult(const Vector &k, Vector &y) const
{
   // Apply the nonlinear form
   pres->SetTime(this->GetTime());
   if (run_mode == 3) {
      tension_func->TryStep(k, dt);
   }
   Hform->Mult(k, y);
}

// Compute the Jacobian from the nonlinear form
Operator &CardiacOperator::GetGradient(const Vector &xp) const
{
   Vector xg;
   /*
   double volume = Hform->GetVolume(xp);
   */
   return Hform->GetGradient(xp);
}

CardiacOperator::~CardiacOperator()
{
   delete J_solver;
   delete fib;
   delete pres;
   if (qat != NULL) {
      delete qat;
   }
   if (tension_func != NULL) {
      delete tension_func;
   }
   delete Q_space;
   if (J_prec != NULL) {
      delete J_prec;
   }
   delete model;
}

// In line visualization
void visualize(ostream &out, ParMesh *mesh, ParGridFunction *deformed_nodes,
               ParGridFunction *field, const char *field_name, bool init_vis)
{
   if (!out)
   {
      return;
   }

   GridFunction *nodes = deformed_nodes;
   int owns_nodes = 0;

   mesh->SwapNodes(nodes, owns_nodes);

   out << "parallel " << mesh->GetNRanks() << " " << mesh->GetMyRank() << "\n";
   out << "solution\n" << *mesh << *field;

   mesh->SwapNodes(nodes, owns_nodes);

   if (init_vis)
   {
      out << "window_size 800 800\n";
      out << "window_title '" << field_name << "'\n";
      if (mesh->SpaceDimension() == 2)
      {
         out << "view 0 0\n"; // view from top
         out << "keys jl\n";  // turn off perspective and light
      }
      out << "keys cm\n";         // show colorbar and mesh
      out << "autoscale value\n"; // update value-range; keep mesh-extents fixed
      out << "pause\n";
   }
   out << flush;
}

