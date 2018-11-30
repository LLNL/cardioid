#include "mfem.hpp"
#include "object.h"
#include "object_cc.hh"
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
#include <dirent.h>
#include <regex.h>
#include <unistd.h>
#include "util.hpp"
#include "MatrixElementPiecewiseCoefficient.hpp"

#define StartTimer(x)
#define EndTimer()

using namespace mfem;

MPI_Comm COMM_LOCAL = MPI_COMM_WORLD;

int main(int argc, char *argv[])
{
   MPI_Init(NULL,NULL);
   int num_ranks, my_rank;
   MPI_Comm_size(COMM_LOCAL,&num_ranks);
   MPI_Comm_rank(COMM_LOCAL,&my_rank);

   std::cout << "Initializing with " << num_ranks << " MPI ranks." << std::endl;

   int order = 1;

   std::vector<std::string> objectFilenames;
   if (argc == 1)
      objectFilenames.push_back("femheard.data");

   for (int iargCursor=1; iargCursor<argc; iargCursor++)
      objectFilenames.push_back(argv[iargCursor]);

   if (my_rank == 0) {
      for (int ii=0; ii<objectFilenames.size(); ii++)
	 object_compilefile(objectFilenames[ii].c_str());
   }
   object_Bcast(0,MPI_COMM_WORLD);

   OBJECT* obj = object_find("femheart", "HEART");
   assert(obj != NULL);

   StartTimer("Read the mesh");
   // Read shared global mesh
   mfem::Mesh *mesh = ecg_readMeshptr(obj, "mesh");
   EndTimer();
   int dim = mesh->Dimension();

   //Fill in the MatrixElementPiecewiseCoefficients
   std::vector<int> heartRegions;
   objectGet(obj,"heart_regions", heartRegions);

   std::vector<double> sigma_i, sigma_e;
   objectGet(obj,"sigma_i",sigma_i);
   objectGet(obj,"sigma_e",sigma_e);

   // Verify Inputs
   assert(heartRegions.size()*3 == sigma_i.size());
   assert(heartRegions.size()*3 == sigma_e.size());

   StartTimer("Constructing the ground part of the mesh.");
   // cardioid_from_ecg = torso/sensor.txt;
   std::unordered_map<int,int> gfFromGid = ecg_readInverseMap(obj,"cardioid_from_ecg");
   EndTimer();
   // Sort+unique pmesh->bdr_attributes and pmesh->attributes?
   StartTimer("Setting Attributes");
   mesh->SetAttributes();
   EndTimer();

   StartTimer("Partition Mesh");
   // If I read correctly, pmeshpart will now point to an integer array
   //  containing a partition ID (rank!) for every element ID.
   int *pmeshpart = mesh->GeneratePartitioning(num_ranks);
   EndTimer();
   
   
   // Elements per rank
   std::map<int,int> local_counts;
   for(int i = 0; i < num_ranks; i++) {
      local_counts[i]=0;
   }

   // Map local indices manually and keep counts for later
   int global_size = mesh->GetNE();
   std::vector<int> local_from_global(global_size);
   std::cout << "Global problem size " << global_size << std::endl;
   for(int i=0; i<global_size; i++) {
      // This calculates everybody's local indices, not just mine.
      local_from_global[i] = local_counts[pmeshpart[i]]++;
   }

   for(int i=0; i<num_ranks; i++) {
      std::cout << "Rank " << i << " has " << local_counts[i] << " elements!" << std::endl;
   }
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh, pmeshpart);

   // Build a new FEC...
   FiniteElementCollection *fec;
   std::cout << "Creating new FEC..." << std::endl;
   fec = new H1_FECollection(order, dim);
   // ...and corresponding FES
   ParFiniteElementSpace *pfespace = new ParFiniteElementSpace(pmesh, fec);
   FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);
   std::cout << "Number of finite element unknowns: "
	     << pfespace->GetTrueVSize() << std::endl;

   // 5. Determine the list of true (i.e. conforming) essential boundary DOFs
   Array<int> ess_tdof_list;   // Essential true degrees of freedom
   // "true" takes into account shared vertices.
   {
      Array<int> ess_bdr(pmesh->bdr_attributes.Max());
      ess_bdr = 0;
      pfespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }

   // 7. Define the solution vector x as a finite element grid function
   //    corresponding to pfespace. Initialize x with initial guess of zero,
   //    which satisfies the boundary conditions.
   ParGridFunction gf_x(pfespace);
   ParGridFunction gf_b(pfespace);
   gf_x = 0.0;
   gf_b = 0.0;

   // Load fiber quaternions from file
   std::shared_ptr<GridFunction> flat_fiber_quat;
   ecg_readGF(obj, "fibers", mesh, flat_fiber_quat);
   std::shared_ptr<ParGridFunction> fiber_quat;
   fiber_quat = std::make_shared<mfem::ParGridFunction>(pmesh, flat_fiber_quat.get(), pmeshpart);

   
   // Load conductivity data?
   MatrixElementPiecewiseCoefficient sigma_m_coeffs(fiber_quat);
   for (int ii=0; ii<heartRegions.size(); ii++) {
      int heartCursor=3*ii;
      Vector sigma_i_vec(&sigma_i[heartCursor],3);
      Vector sigma_e_vec(&sigma_e[heartCursor],3);
      Vector sigma_m_vec(3);
      for (int jj=0; jj<3; jj++)
      {
         sigma_m_vec[jj] = (sigma_e_vec[jj]*sigma_i_vec[jj])/(sigma_i_vec[jj]+sigma_e_vec[jj]);
      }
    
      sigma_m_coeffs.heartConductivities_[heartRegions[ii]] = sigma_m_vec;
   }

   // 8. Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the Laplacian operator -Delta, by adding the Diffusion
   //    domain integrator.

   StartTimer("Forming bilinear system (heart)");

   ParBilinearForm *b = new ParBilinearForm(pfespace);
   b->AddDomainIntegrator(new DiffusionIntegrator(sigma_m_coeffs));
   b->Assemble();
   // This creates the linear algebra problem.
   HypreParMatrix heart_mat;
   b->FormSystemMatrix(ess_tdof_list, heart_mat);
   // Parallel versions don't get finalized? // heart_mat.Finalize();

   EndTimer();

   StartTimer("Forming bilinear system (torso)");
   
   // Brought out of loop to avoid unnecessary duplication
   ParBilinearForm *a = new ParBilinearForm(pfespace);   // defines a.
   // this is the Laplacian: grad u . grad v with linear coefficient.
   // we defined "one" ourselves in step 6.
   a->AddDomainIntegrator(new DiffusionIntegrator(sigma_m_coeffs));
   // a->Assemble();   // This creates the loops.

   a->Update(pfespace);
   a->Assemble();
   HypreParMatrix torso_mat;
   a->FormSystemMatrix(ess_tdof_list,torso_mat);
   // Look into FLS to see what steps are being applied and make sure boundary conditions are still being set
   
   // 10. Define a simple symmetric Gauss-Seidel preconditioner and use it to
   //     solve the system A X = B with PCG.
   EndTimer();

   //FIXME!!!  move me outside the loop
   HyprePCG pcg(torso_mat);
   pcg.SetTol(1e-12);
   pcg.SetMaxIter(2000);
   pcg.SetPrintLevel(2);
   HypreSolver *M_test = new HypreBoomerAMG(torso_mat);
   pcg.SetPreconditioner(*M_test);
   //end move me

   Vector phi_e;
   Vector phi_b;
   

#if 0
      //Do we have a Vm file present?
      if (access((VmFilename + "000000").c_str(), R_OK) == -1) { continue; }

      std::shared_ptr<ParGridFunction> gf_Vm;
      std::shared_ptr<GridFunction> flat_gf_Vm;
      std::vector<double> gfvmData;
      double time = ecg_readParGF(obj, VmFilename, global_size, gfFromGid, gfvmData);
      flat_gf_Vm = std::make_shared<mfem::GridFunction>(fespace, gfvmData.data());
      gf_Vm = std::make_shared<mfem::ParGridFunction>(pmesh, flat_gf_Vm.get(), pmeshpart);

      StartTimer("Solve");
      
      heart_mat.Mult(*gf_Vm, gf_b);
      a->FormLinearSystem(ess_tdof_list,gf_x,gf_b,torso_mat,phi_e,phi_b);


      //HypreSmoother M(torso_mat, 6 /*GS*/);
      // PCG(torso_mat, M, phi_b, phi_e, 1, 2000, 1e-12);//, 0.0);
      phi_b *= -1.0;
      pcg.Mult(phi_b,phi_e);

      EndTimer();
      
      // 11. Recover the solution as a finite element grid function.
      a->RecoverFEMSolution(phi_e, phi_b, gf_x);

      // 12. Save the refined mesh and the solution. This output can be viewed later
      //     using GLVis: "glvis -m refined.mesh -g sol.gf".
      for (int ielec=0; ielec<fileFromElectrode.size(); ielec++)
      {
	 if(my_rank == pmeshpart[gfidFromElectrode[ielec]]) {
            //CHECKME, can I do this access?!
            fileFromElectrode[ielec] << time << "\t" << gf_x[local_from_global[gfidFromElectrode[ielec]]] << std::endl;
         }
      }
      
#ifdef DEBUG
      std::ofstream sol_ofs(outDir+"/sol"+std::to_string(time)+".gf");
      sol_ofs.precision(8);
      gf_x.SaveAsOne(sol_ofs);
      sol_ofs.close();
#endif	 
   }

#endif

#ifdef DEBUG
   std::ofstream mesh_ofs(outDir+"/refined.mesh");
   mesh_ofs.precision(8);
   pmesh->PrintAsOne(mesh_ofs);
#endif
   
   // 14. Free the used memory.
   delete M_test;
   delete a;
   delete b;
   delete pfespace;
   if (order > 0) { delete fec; }
   delete mesh, pmesh, pmeshpart;
   
   return 0;
}
