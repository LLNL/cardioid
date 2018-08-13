#include "mfem.hpp"
#include "object.h"
#include "ddcMalloc.h"
#include "pio.h"
#include "pioFixedRecordHelper.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <cassert>
#include <memory>
#include <set>
#include <ctime>
#include "ecgdefs.hpp"
#include "ecgobjutil.hpp"

using namespace mfem;

MPI_Comm COMM_LOCAL = MPI_COMM_WORLD;

int main(int argc, char *argv[])
{
  MPI_Init(NULL,NULL);
  // 1. Parse command-line options.
  int order = 1;
  time_t timestamp;

  ecg_process_args(argc,argv);

  /* Inputs:
     - "mesh" as .vtk (example: slab.vtk)
     - points as float[9] (143871)
     - cells as int[9] (135000x9=1215000)
     - cell types as int (135000)
     - materials as int[9] (135000)
     - "fibers" as .gf (example: slab.FiberQuat.gf)
     NOTE: header also contains FEC type, L2_3D_P0
     - values as int[3] (135000)
     - "Vm" as .gf (example: slab.Vm.gf)
     NOTE: header also contains FEC type, H1_3D_P1
     - values as int (63240)
     - "ground" as .set (example: slab.ground.set)
     - Node IDs as int (1581)
  */

  OBJECT* obj = object_find("ecg", "ECG");
  assert(obj != NULL);

  // Read unordered set of grounds
  std::set<int> ground = ecg_readSet(obj, "ground");

  // Read shared global mesh
  Mesh *mesh = ecg_readMeshptr(obj, "mesh");
  int dim = mesh->Dimension();

  // KNOWNS thus far:  mesh, dim, ground

  //Fill in the MatrixElementPiecewiseCoefficients
  std::vector<int> bathRegions;
  objectGetv(obj,"bath_regions",bathRegions);
  std::vector<int> heartRegions;
  objectGetv(obj,"heart_regions", heartRegions);
  std::vector<double> sigma_b;
  objectGetv(obj,"sigma_b",sigma_b);;
  std::vector<double> sigma_i;
  std::vector<double> sigma_e;
  objectGetv(obj,"sigma_i",sigma_i);
  objectGetv(obj,"sigma_e",sigma_e);

  // Verify Inputs
  assert(bathRegions.size() == sigma_b.size());
  assert(heartRegions.size()*3 == sigma_i.size());
  assert(heartRegions.size()*3 == sigma_e.size());

  // Added KNOWNS: bathRegions, heartRegions, sigma_{b,i,e}

  // Read sensors
  std::unordered_map<int,int> sensors = ecg_readInverseMap(obj,"sensors");

  // Read gid->gf mapping

  // Describe input data
#ifdef DEBUG
  std::cout << "Read " << mesh->GetNE() << " elements, "
	    << mesh->GetNV() << " vertices, "
	    << mesh->GetNBE() << " boundary elements, and "
	    << mesh->GetNumFaces() << " faces in "
	    << dim << " dimensions." << std::endl;
#endif

  {
    // Iterate over boundary elements
    int nbe=mesh->GetNBE();
    for(int i=0; i<nbe; i++){
      Element *ele = mesh->GetBdrElement(i);
      const int *v = ele->GetVertices();
      const int nv = ele->GetNVertices();
      // Search for element's vertices in the ground set
      bool isGround = true;
      for( int ivert=0; ivert<nv; ivert++) {
	if (ground.find(v[ivert]) == ground.end()) {
	  isGround = false;
	  break;
	}
      }
      // Set region type accordingly per ecg.data
      if (isGround) {
	ele->SetAttribute(2); // ess[1] below
      } else {
	ele->SetAttribute(1); // ess[0] below
      }
    }
  }
  // Sort+unique mesh->bdr_attributes and mesh->attributes?
  mesh->SetAttributes();

  // Pointer to the FEC stored in mesh's internal GF
  FiniteElementCollection *fec;
  if (mesh->GetNodes()) {
    // Not called?
    fec = mesh->GetNodes()->OwnFEC();
    std::cout << "Using isoparametric FEs: " << fec->Name() << std::endl;
  }
  else {
    std::cout << "Creating new FEC..." << std::endl;
    fec = new H1_FECollection(order, dim);
  }
  
  FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);
  std::cout << "Number of finite element unknowns: "
	    << fespace->GetTrueVSize() << std::endl;

  // GetTrueVSize() gets all vector-local dofs
  // GetGlobalTDofNumber converts to a GF index reference against gidfromgf
   
  // 5. Determine the list of true (i.e. conforming) essential boundary
  //    dofs.
  Array<int> ess_tdof_list;   // Essential true degrees of freedom
  // "true" takes into account shared vertices.
  if (mesh->bdr_attributes.Size()) {
    assert(mesh->bdr_attributes.Max() > 1 && "Can't find a ground boundary!");
    Array<int> ess_bdr(mesh->bdr_attributes.Max());
    ess_bdr = 0;
    ess_bdr[1] = 1;
    fespace->GetEssentialTrueDofs(ess_bdr /*const*/, ess_tdof_list);
  }
   

  // 7. Define the solution vector x as a finite element grid function
  //    corresponding to fespace. Initialize x with initial guess of zero,
  //    which satisfies the boundary conditions.
  GridFunction gf_x(fespace);
  GridFunction gf_b(fespace);
  gf_x = 0.0;
   

  // Load fiber quaterions from file
  std::shared_ptr<GridFunction> fiber_quat;
  ecg_readGF(obj, "fibers", mesh, fiber_quat);

  // Load conductivity data?
  MatrixElementPiecewiseCoefficient sigma_i_coeffs(fiber_quat);
  MatrixElementPiecewiseCoefficient sigma_ie_coeffs(fiber_quat);
  for (int ii=0; ii<heartRegions.size(); ii++) {
    int heartCursor=3*ii;
    Vector sigma_i_vec(&sigma_i[heartCursor],3);
    Vector sigma_e_vec(&sigma_e[heartCursor],3);
    Vector sigma_ie_vec = sigma_e_vec;
    sigma_ie_vec += sigma_i_vec;
    
    sigma_i_coeffs.heartConductivities_[heartRegions[ii]] = sigma_i_vec;
    sigma_ie_coeffs.heartConductivities_[heartRegions[ii]] = sigma_ie_vec;
  }
  for (int ii=0; ii<bathRegions.size(); ii++) {
    sigma_i_coeffs.bathConductivities_[bathRegions[ii]] = 0;
    sigma_ie_coeffs.bathConductivities_[bathRegions[ii]] = sigma_b[ii];
  }

  // 8. Set up the bilinear form a(.,.) on the finite element space
  //    corresponding to the Laplacian operator -Delta, by adding the Diffusion
  //    domain integrator.

#ifdef DEBUG
  timestamp = time(nullptr);
  std::cout << std::endl << "Forming bilinear system (heart)... ";
#endif
  BilinearForm *b = new BilinearForm(fespace);
  b->AddDomainIntegrator(new DiffusionIntegrator(sigma_i_coeffs));
  b->Assemble();
  // This creates the linear algebra problem.
  SparseMatrix heart_mat;
  b->FormSystemMatrix(ess_tdof_list, heart_mat);
  heart_mat.Finalize();
#ifdef DEBUG
  std::cout << "Done in " << time(nullptr)-timestamp << "s." << std::endl;
#endif

  // cardioid_from_ecg = torso/sensor.txt;
  std::unordered_map<int,int> gfFromGid = ecg_readInverseMap(obj,"cardioid_from_ecg");

  // Want to use the below instead of this
  std::string VmFilename;
  objectGet(obj, "Vm", VmFilename, ""); // Vm = ../test.Vm.gf;

  // Want to use this instead of the above
  std::string VmPattern;
  objectGet(obj, "VmPattern", VmPattern, ""); // VmPattern = ../torsoRun/snapshot.%012d/Vm#%06d;
  
  std::shared_ptr<GridFunction> gf_Vm;
  {
    //char *VmFileCstr = new char[VmPattern.length()+1]; // Cannot use variable formats!
    //sprintf(VmFileCstr, VmPattern.c_str(), 200, i);
    PFILE* file = Popen(VmFilename.c_str(), "r", COMM_LOCAL);
    gf_Vm = std::make_shared<mfem::GridFunction>(fespace);
    
    OBJECT* hObj = file->headerObject;
    std::vector<std::string> fieldNames,  fieldTypes;
    objectGetv(hObj, "field_names", fieldNames); // field_names = gid Vm dVmD dVmR;
    objectGetv(hObj, "field_types", fieldTypes); // field_types = u f f f;
    assert(file->datatype == FIXRECORDASCII);
    
    PIO_FIXED_RECORD_HELPER* helper = (PIO_FIXED_RECORD_HELPER*) file->helper;
    unsigned lrec = helper->lrec;
    unsigned nRecords = file->bufsize/lrec;
    for (unsigned int irec=0; irec<nRecords; irec++) { //read in the file
      unsigned maxRec = 2048;          // Maximum record *length*
      char buf[maxRec+1];
	
      Pfgets(buf, maxRec, file);
      assert(strlen(buf) < maxRec);

      std::stringstream readline(buf);

      int gid;   readline >> gid; // For now we just happen to know gid is column 1
      double Vm; readline >> Vm;  // For now we just happen to know the data is column 2

      (*gf_Vm)[gfFromGid[gid]] = Vm;
    }
      
    Pclose(file);
    //}
  }

  heart_mat.Mult(*gf_Vm, gf_b);
#ifdef DEBUG
  std::cout << "heart_mat.Empty()" << " "
	    << "heart_mat.Finalized()" << " "
	    << "heart_mat.MaxNorm()" << " "
	    << "heart_mat.Height()" << " "
	    << "heart_mat.Width()" << " "
	    << "heart_mat.ActualWidth()" << std::endl;
  std::cout << heart_mat.Empty() << " "
	    << heart_mat.Finalized() << " "
	    << heart_mat.MaxNorm() << " "
	    << heart_mat.Height() << " "
	    << heart_mat.Width() << " "
	    << heart_mat.ActualWidth() << std::endl;
  std::cout << "gf_Vm->Max()" << " " << "gf_Vm->Min()" << std::endl;
  std::cout << gf_Vm->Max() << " " << gf_Vm->Min() << std::endl;
  std::cout << "gf_b.Max()" << " " << "gf_b.Min()" << std::endl;
  std::cout << gf_b.Max() << " " << gf_b.Min() << std::endl;
#endif

#ifdef DEBUG
  timestamp = time(nullptr);
  std::cout << std::endl << "Forming bilinear system (torso)... ";
#endif
  BilinearForm *a = new BilinearForm(fespace);   // defines a.
  // this is the Laplacian: grad u . grad v with linear coefficient.
  // we defined "one" ourselves in step 6.
  a->AddDomainIntegrator(new DiffusionIntegrator(sigma_ie_coeffs));
  a->Assemble();   // This creates the loops.
  SparseMatrix torso_mat;
  Vector phi_e;
  Vector phi_b;
  a->FormLinearSystem(ess_tdof_list,gf_x,gf_b,torso_mat,phi_e,phi_b);
#ifdef DEBUG
  std::cout << "Done in " << time(nullptr)-timestamp << "s." << std::endl;
#endif
   
  // NOTE THE ifdef
#ifndef MFEM_USE_SUITESPARSE
  // 10. Define a simple symmetric Gauss-Seidel preconditioner and use it to
  //     solve the system A X = B with PCG.
#ifdef DEBUG
  timestamp = time(nullptr);
  std::cout << std::endl << "Solving... ";
#endif
  GSSmoother M(torso_mat);
  PCG(torso_mat, M, phi_b, phi_e, 1, 2000, 1e-12, 0.0);
#ifdef DEBUG
  std::cout << "Done in " << time(nullptr)-timestamp << "s." << std::endl;
#endif
#else
  // 10. If MFEM was compiled with SuiteSparse, use UMFPACK to solve the system.
  UMFPackSolver umf_solver;
  umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
  umf_solver.SetOperator(torso_mat);
  umf_solver.Mult(phi_b, phi_e);
  // See parallel version for HypreSolver - which is an LLNL package.
#endif

#ifdef DEBUG
  std::cout << std::endl << "phi_e.Max()" << " "
	    << "phi_e.Min()" << " "
	    << "phi_e.Norml2()" << " "
	    << "phi_e.Size():" << std::endl
	    << phi_e.Max() << " "
	    << phi_e.Min() << " "
	    << phi_e.Norml2() << " "
	    << phi_e.Size() << std::endl;
#endif

  // 11. Recover the solution as a finite element grid function.
  a->RecoverFEMSolution(phi_e, phi_b, gf_x);

#ifdef DEBUG
  std::cout << std::endl << "gf_x.Max()" << " "
	    << "gf_x.Min()" << " "
	    << "gf_x.Norml2()" << " "
	    << "gf_x.Size()" << std::endl;
  std::cout << gf_x.Max() << " "
	    << gf_x.Min() << " "
	    << gf_x.Norml2() << " "
	    << gf_x.Size() << std::endl;
#endif
   
  // 12. Save the refined mesh and the solution. This output can be viewed later
  //     using GLVis: "glvis -m refined.mesh -g sol.gf".
  std::ofstream mesh_ofs("refined.mesh");
  mesh_ofs.precision(8);
  mesh->Print(mesh_ofs);
  std::ofstream sol_ofs("sol.gf");
  std::ofstream sol_vtk("sol.vtk");
  sol_ofs.precision(8);
  gf_x.Save(sol_ofs);
  mesh->PrintVTK(sol_vtk);
  gf_x.SaveVTK(sol_vtk,"test",1);

  // 14. Free the used memory.
  delete a;
  delete b;
  delete fespace;
  if (order > 0) { delete fec; }
  delete mesh;

  return 0;
}
