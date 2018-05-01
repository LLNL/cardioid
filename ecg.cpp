
#include "mfem.hpp"
#include "object.h"
#include "ddcMalloc.h"
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <cassert>

using namespace std;
using namespace mfem;

MPI_Comm COMM_LOCAL = MPI_COMM_WORLD;

DenseMatrix quat2rot(const Vector& quat);

class MatrixElementPiecewiseCoefficient : public MatrixCoefficient
{
 public:
   MatrixElementPiecewiseCoefficient() : MatrixCoefficient(3) {}

   virtual void Eval(DenseMatrix &K, ElementTransformation& T, const IntegrationPoint &ip)
   {
      unordered_map<int,Vector>::iterator iter = heartConductivities_.find(T.Attribute);
      if (iter != heartConductivities_.end())
      {
         Vector direction(3);
         p_gf_->GetVectorValue(T.ElementNo, ip, direction);
         Vector quat(4);
         double w2 = 1;
         for (int ii=0; ii<3; ii++)
         {
            quat(ii+1) = direction(ii);
            w2 -= direction(ii)*direction(ii);
         }
         quat(0) = sqrt(w2);
         
         DenseMatrix VVV = quat2rot(quat);
         MultADAt(VVV,iter->second,K);
      }
      else
      {
         K=0.0;
         unordered_map<int,double>::iterator iter = bathConductivities_.find(T.Attribute);
         if (iter != bathConductivities_.end())
         {
            for (int ii=0; ii<3; ii++)
            {
               K(ii,ii) = iter->second;
            }
         }
      }
   }

   shared_ptr<GridFunction> p_gf_;
   unordered_map<int,Vector> heartConductivities_;
   unordered_map<int,double> bathConductivities_;
};


/** The object_get code doesn't handle empty strings for a default value
 *  very well.  The problem is that object_parse attempts to tokenize
 *  the default value and strtok_r returns a NULL pointer.  This results
 *  in a situation where the pointer passed to object_get (i.e., tmp)
 *  isn't set by object_get.  Fortunately, this can be detected by the
 *  fact that object_get will return 0 instead of 1.  We can catch this
 *  case and do the right thing with the value we return to the caller.
 */
void objectGet(OBJECT* obj, const string& name, string& value, const string& defVal)
{
   char* tmp;
   int nFound = object_get(obj, name.c_str(), &tmp, STRING, 1, defVal.c_str());
   if (nFound != 0)
   {
      value = tmp;
      ddcFree(tmp);
   }
   else
      value = "";
}


int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   int order = 1;

   vector<string> objectFilenames;
   if (argc == 1)
   {
      objectFilenames.push_back("ecg.data");
   }
   for (int iargCursor=1; iargCursor<argc; iargCursor++)
   {
      objectFilenames.push_back(argv[iargCursor]);
   }
   for (int ii=0; ii<objectFilenames.size(); ii++)
   {
      object_compilefile(objectFilenames[ii].c_str());
   }

   OBJECT* obj = object_find("ecg", "ECG");
   assert(obj != NULL);
   string mesh_file;
   objectGet(obj, "mesh", mesh_file, "");

   Mesh *mesh = new Mesh(mesh_file.c_str(), 1, 1);
   int dim = mesh->Dimension();

   FiniteElementCollection *fec;
   if (mesh->GetNodes())
   {
      fec = mesh->GetNodes()->OwnFEC();
      cout << "Using isoparametric FEs: " << fec->Name() << endl;
   }
   else
   {
      fec = new H1_FECollection(order, dim);
   }
   FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);
   cout << "Number of finite element unknowns: "
        << fespace->GetTrueVSize() << endl;

   // 5. Determine the list of true (i.e. conforming) essential boundary
   //    dofs.
   //    In this example, the boundary conditions are defined by marking all
   //    the boundary attributes from the mesh as essential (Dirichlet) and
   //    converting them to a list of true dofs.
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
   LinearForm *b = new LinearForm(fespace);
   ConstantCoefficient one(1.0);  // coef in front of grad u grad v
                                  // We'll change to tensor coefficient
                                  // (sigma -- Jamie will look for this in
                                  // MFEM, or add.
   b->AddDomainIntegrator(new DomainLFIntegrator(one));
   b->Assemble();

   // 7. Define the solution vector x as a finite element grid function
   //    corresponding to fespace. Initialize x with initial guess of zero,
   //    which satisfies the boundary conditions.
   GridFunction x(fespace);
   x = 0.0;   // essential boundary conditions are zero, so set whole thing
              // to zero.
   

   //Fill in the MatrixElementPiecewiseCoefficients
   MatrixElementPiecewiseCoefficient sigma_i_coeffs;
   MatrixElementPiecewiseCoefficient sigma_ie_coeffs;

   // 8. Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the Laplacian operator -Delta, by adding the Diffusion
   //    domain integrator.
   BilinearForm *a = new BilinearForm(fespace);   // defines a.
   // this is the Laplacian: grad u . grad v with linear coefficient.
   // we defined "one" ourselves in step 6.
   a->AddDomainIntegrator(new DiffusionIntegrator(sigma_ie_coeffs));
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

   // 12. Save the refined mesh and the solution. This output can be viewed later
   //     using GLVis: "glvis -m refined.mesh -g sol.gf".
   ofstream mesh_ofs("refined.mesh");
   mesh_ofs.precision(8);
   mesh->Print(mesh_ofs);
   ofstream sol_ofs("sol.gf");
   sol_ofs.precision(8);
   x.Save(sol_ofs);

   // 14. Free the used memory.
   delete a;
   delete b;
   delete fespace;
   if (order > 0) { delete fec; }
   delete mesh;

   return 0;
}
