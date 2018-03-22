// ecg - heart and torso
//
// Compile with: make ecg
//
// Run:          ecg -m <mesh_file>

// Modified from mfem/examples/ex1.cpp.

#include "mfem.hpp"
#include "pio.h"

#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

MPI_Comm COMM_LOCAL = MPI_COMM_WORLD;

void sigma_fct (const Vector &x, DenseMatrix &f);

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = "../data/star.mesh";
   int order = 1;
   bool static_cond = false;
   bool visualization = 1;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree) or -1 for"
                  " isoparametric space.");
   args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
                  "--no-static-condensation", "Enable static condensation.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   // 2. Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral, hexahedral, surface and volume meshes with
   //    the same code.
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();

   // DKTMP
   int n_attributes = mesh->attributes.Max ();
   cerr << "n_attributes: " << n_attributes << endl;

   // 3. Refine the mesh to increase the resolution. In this example we do
   //    'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
   //    largest number that gives a final mesh with no more than 50,000
   //    elements.
   {
      int ref_levels =
         (int)floor(log(50000./mesh->GetNE())/log(2.)/dim);
      for (int l = 0; l < ref_levels; l++)
      {
         mesh->UniformRefinement();
      }
   }

   // 4. Define a finite element space on the mesh. Here we use continuous
   //    Lagrange finite elements of the specified order. If order < 1, we
   //    instead use an isoparametric/isogeometric space.
   //    dk: order = no. polynomial basis functions per element (scaled by
   //                           dimension)
   //                e.g., 1D order 1 = 2 linear functions per element
   //                      if 1D with a node in between: still order 1.
   //        defining master element (unit sides, etc.)
   FiniteElementCollection *fec;
   if (order > 0)
   {
      fec = new H1_FECollection(order, dim);
      // in our case: order = 1, dim = 3.  order is the degree of the basis
      //                                   fct.
   }
   else if (mesh->GetNodes())
   {
      fec = mesh->GetNodes()->OwnFEC();
      cout << "Using isoparametric FEs: " << fec->Name() << endl;
   }
   else
   {
      fec = new H1_FECollection(order = 1, dim);
   }
   // dk: vtk file, for example -- discretizes,  Also give it "master
   //     element" template.
   // It comes up with no. unknowns, mapping between physical elements and
   // master (how transformed the unit-side element is).
   // one tet = 4 unknowns (one per vertex).  Two tetrahedrons face-to-face
   // = 5.
   FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);
   cout << "Number of finite element unknowns: "
        << fespace->GetTrueVSize() << endl;

   // 5. Determine the list of true (i.e. conforming) essential boundary
   //    dofs.
   //    In this example, the boundary conditions are defined by marking all
   //    the boundary attributes from the mesh as essential (Dirichlet) and
   //    converting them to a list of true dofs.
   // dk: mark faces.  Dirichlet = known concentration.  A Neumann flux
   // condition is not setting a value; instead, add a term to system
   // (later).
   Array<int> ess_tdof_list;   // Essential true degrees of freedom
                               // "true" takes into account shared vertices.
   if (mesh->bdr_attributes.Size())
   {
      Array<int> ess_bdr(mesh->bdr_attributes.Max());
      ess_bdr = 1;
      fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }

   // 6. Set up the linear form b(.) which corresponds to the right-hand side of
   //    the FEM linear system, which in this case is (1,phi_i) where phi_i are
   //    the basis functions in the finite element fespace.
   // dk: if had a forcing function (source term) on RHS.
   //     In our case there is none.  So replace these (later) in linear object.
   LinearForm *b = new LinearForm(fespace);

   // Note: this defines variable "one".
   ConstantCoefficient one(1.0);

   b->AddDomainIntegrator(new DomainLFIntegrator(one));
   b->Assemble();

   // 7. Define the solution vector x as a finite element grid function
   //    corresponding to fespace. Initialize x with initial guess of zero,
   //    which satisfies the boundary conditions.
   // dk: grid fct is vector with extra info.  vec of length of no. unknowns
   // in system -- extra info: map to physical space; connectivity/shared
   // vertices; knows how to interpolate.
   GridFunction x(fespace);
   x = 0.0;   // essential boundary conditions are zero, so set whole thing
              // to zero.

   // 8. Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the Laplacian operator -Delta, by adding the Diffusion
   //    domain integrator.
   // dk: setting up LHS.
   //     "test function" "trial function" -- multiply A matrix by test and
   //     trial function.  "Bi-linear" - linear in each of the terms.
   //     "v" and "u".  (Not really "test" and "trial" because linear algebra
   //     does it.)
   //     a is finite-element-object with instructions on how to form matrix.
   BilinearForm *a = new BilinearForm(fespace);   // defines a.

   // _________________________________
   // Conductance: scalar, depending on tissue type outside heart; tensor within
   // heart.  Integrator for each, added; must return zeros outside their
   // respective domains.

   // .................................
   // Tensor conductance.

   /* Constant-matrix test.

   // DenseMatrix (n): square matrix of size n.
   DenseMatrix sigma_3by3(3);

   // 3 x 3 matrix, given by columns.
   double elem_array[] = {1.0, 0.0, 0.0,
                          0.0, 1.0, 0.0,
                          0.0, 0.0, 1.0};
   sigma_3by3 = elem_array;

   // DKTMP
   cerr << "sigma_3by3[0][1]: " << sigma_3by3.Elem (0, 1) << endl;
   cerr << "sigma_3by3: " << endl;
   sigma_3by3.Print (cerr);
   */

   // MatrixFunctionCoefficient ().  See fem/coefficient.hpp line 476ff.
   // sigma_fct () defined below.
   MatrixFunctionCoefficient sigma (dim, sigma_fct);

   //DKTMP Alt:
   //MatrixConstantCoefficient sigma (sigma_3by3);

   // .................................
   // Scalar conductance.
   int n_tissue_types = mesh->attributes.Max ();
   Vector k(n_tissue_types);

   // DKTMP: need to set k_vals(0), k(1)...
   k = 0.0;
   k(1) = 5.0;    // 3 x 3 x 3 - hex.small.vtk -- center cell.
   k(0) = 0.0;
   PWConstCoefficient k_fct (k);

   // .................................
   // Add both integrators.
   a->AddDomainIntegrator(new DiffusionIntegrator (sigma));
   a->AddDomainIntegrator(new DiffusionIntegrator (k_fct));
   // _________________________________

   // 9. Assemble the bilinear form and the corresponding linear system,
   //    applying any necessary transformations such as: eliminating boundary
   //    conditions, applying conforming constraints for non-conforming AMR,
   //    static condensation, etc.
   // dk: at this point, get out linear object (linear system)
   //     this sets up the loops Rob talked about. When we set up the
   //     domain integrator, above, it now knows to insert that in the
   //     innermost loop.
   if (static_cond) { a->EnableStaticCondensation(); }
   a->Assemble();   // This creates the loops.

   SparseMatrix A;   // This is what we want.
   Vector B, X;
   // This creates the linear algebra problem.
   a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);

   cout << "Size of linear system: " << A.Height() << endl;
   // true dof minus essential unknowns (we defined as known).

// NOTE THE ifdef
#ifndef MFEM_USE_SUITESPARSE
   // 10. Define a simple symmetric Gauss-Seidel preconditioner and use it to
   //     solve the system A X = B with PCG.
   // dk: can use either preconditioned conjugate gradient (simplest linear
   //     solver); or UMPFACK - direct solver (gaussian elimination).
   //     (PCG better for giant systems, sparse matrices.  HPC > 99% of
   //     time use iterative method.)
   GSSmoother M(A);
   PCG(A, M, B, X, 1, 200, 1e-12, 0.0);
#else
   // 10. If MFEM was compiled with SuiteSparse, use UMFPACK to solve the system.
   UMFPackSolver umf_solver;
   umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
   umf_solver.SetOperator(A);
   umf_solver.Mult(B, X);
   // See parallel version for HypreSolver - which is an LLNL package.
#endif

   // 11. Recover the solution as a finite element grid function.
   a->RecoverFEMSolution(X, *b, x);

   // 12. Save the refined mesh and the solution. This output can be viewed
   //     later using GLVis: "glvis -m refined.mesh -g sol.gf".
   ofstream mesh_ofs("refined.mesh");
   mesh_ofs.precision(8);
   mesh->Print(mesh_ofs);
   ofstream sol_ofs("sol.gf");
   sol_ofs.precision(8);
   x.Save(sol_ofs);

   // 13. Send the solution by socket to a GLVis server.
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock.precision(8);
      sol_sock << "solution\n" << *mesh << x << flush;
   }

   // 14. Free the used memory.
   delete a;
   delete b;
   delete fespace;
   if (order > 0) { delete fec; }
   delete mesh;

   return 0;
}

int sigma_calls = 0;
void sigma_fct (const Vector &x, DenseMatrix &f)
{
   //DKTMP
   //cerr << "[sigma_fct] x: " << x(0) << " " << x(1) << " " << x(2) << endl;

   // DenseMatrix (n): square matrix of size n.
   DenseMatrix sigma_3by3(3);

   // Create diagonal matrix (3 x 3, with diagonal elements = 1.0)
   // 3 x 3 matrix, given by columns.
   double elem_array[] = {1.00, 0.0,  0.0,
                          0.0,  1.00, 0.0,
                          0.0,  0.0,  0.01};
   sigma_3by3 = elem_array;

   f = sigma_3by3;

   /*
   sigma_calls++;
   if (sigma_calls > 5) {
      cerr << "Exiting in sigma_fct" << endl;
      exit (1);
   }
   */
}
