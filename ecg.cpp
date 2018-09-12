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
#include <dirent.h>
#include <regex.h>
#include "ecgdefs.hpp"
#include "ecgobjutil.hpp"

using namespace mfem;

MPI_Comm COMM_LOCAL = MPI_COMM_WORLD;

int main(int argc, char *argv[])
{
   MPI_Init(NULL,NULL);
   int num_ranks;
   MPI_Comm_size(COMM_LOCAL,&num_ranks);

   std::cout << "Initializing with " << num_ranks << " MPI ranks." << std::endl;
   // 1. Parse command-line options.
   int order = 1;

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

   StartTimer("Read the mesh");
   // Read shared global mesh
   mfem::Mesh *mesh = ecg_readMeshptr(obj, "mesh");
   EndTimer();
   int dim = mesh->Dimension();
   StartTimer("Partition Mesh");
   int *pmeshpart = mesh->GeneratePartitioning(num_ranks);
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);//, pmeshpart);
   EndTimer();
   
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

   StartTimer("Constructing the ground part of the mesh.");
   // cardioid_from_ecg = torso/sensor.txt;
   std::unordered_map<int,int> gfFromGid = ecg_readInverseMap(obj,"cardioid_from_ecg");

   {
      // Iterate over boundary elements
      int nbe=pmesh->GetNBE();
      for(int i=0; i<nbe; i++){
	 Element *ele = pmesh->GetBdrElement(i);
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
   EndTimer();
   // Sort+unique pmesh->bdr_attributes and pmesh->attributes?
   StartTimer("Setting Attributes");
   pmesh->SetAttributes();
   EndTimer();
   
   // Pointer to the FEC stored in mesh's internal GF
   FiniteElementCollection *fec;
   if (pmesh->GetNodes()) {
      // Not called?
      fec = pmesh->GetNodes()->OwnFEC();
      std::cout << "Using isoparametric FEs: " << fec->Name() << std::endl;
   }
   else {
      std::cout << "Creating new FEC..." << std::endl;
      fec = new H1_FECollection(order, dim);
   }
  
   ParFiniteElementSpace *pfespace = new ParFiniteElementSpace(pmesh, fec);
   std::cout << "Number of finite element unknowns: "
	     << pfespace->GetTrueVSize() << std::endl;

   // GetTrueVSize() gets all vector-local dofs
   // GetGlobalTDofNumber converts to a GF index reference against gidfromgf
   
   // 5. Determine the list of true (i.e. conforming) essential boundary
   //    dofs.
   Array<int> ess_tdof_list;   // Essential true degrees of freedom
   // "true" takes into account shared vertices.
   if (pmesh->bdr_attributes.Size()) {
      assert(pmesh->bdr_attributes.Max() > 1 && "Can't find a ground boundary!");
      Array<int> ess_bdr(pmesh->bdr_attributes.Max());
      ess_bdr = 0;
      ess_bdr[1] = 1;
      pfespace->GetEssentialTrueDofs(ess_bdr /*const*/, ess_tdof_list);
   }


   // 7. Define the solution vector x as a finite element grid function
   //    corresponding to pfespace. Initialize x with initial guess of zero,
   //    which satisfies the boundary conditions.
   ParGridFunction gf_x(pfespace);
   ParGridFunction gf_b(pfespace);
   gf_x = 0.0;


   // Load fiber quaternions from file
   std::shared_ptr<ParGridFunction> fiber_quat;
   ecg_readGF(obj, "fibers", pmesh, fiber_quat);

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

   StartTimer("Forming bilinear system (heart)");

   ParBilinearForm *b = new ParBilinearForm(pfespace);
   b->AddDomainIntegrator(new DiffusionIntegrator(sigma_i_coeffs));
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
   a->AddDomainIntegrator(new DiffusionIntegrator(sigma_ie_coeffs));
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
   
   // Read in the electrode list
   std::string electrodeFilename;
   objectGet(obj, "ecg_electrodes", electrodeFilename, "");
   std::vector<std::string> nameFromElectrode;
   std::vector<int> gfidFromElectrode;
   {
      std::ifstream electrodeFile(electrodeFilename.c_str());
      while (!!electrodeFile)
      {
         std::string name;
         std::string dontCare;
         int gfid;
         electrodeFile >> name;
         electrodeFile >> dontCare;
         electrodeFile >> gfid;

         if (!electrodeFile) { break; }
         nameFromElectrode.push_back(name);
         gfidFromElectrode.push_back(gfid);
      }
   }
   // Want to use this instead of literal filenames for multiple time steps
   std::string VmSubfile;
   objectGet(obj, "vm_subfile", VmSubfile, ""); // VmPattern = ../torsoRun/snapshot.%012d/Vm#%06d;
   std::string rootFilename;
   objectGet(obj, "simdir", rootFilename, ".");
   std::string outDir;
   objectGet(obj, "outdir", outDir, rootFilename.c_str());
   std::vector<std::ofstream> fileFromElectrode(nameFromElectrode.size());
   for (int ielec=0; ielec<fileFromElectrode.size(); ielec++)
   {
      fileFromElectrode[ielec] = std::ofstream(outDir+"/"+nameFromElectrode[ielec]+".txt");
   }
   
   DIR *dir;
   dir = opendir(rootFilename.c_str());
   if (dir == NULL) return 0;

   dirent *entry;
   regex_t snapshotRegex;
   int retCode = regcomp(&snapshotRegex, "^snapshot\\.[[:digit:]]\\{1,\\}$", REG_NOSUB);
   assert(retCode == 0);
   while((entry = readdir(dir)) != NULL)
   {
      //Does the file match the output pattern?
      retCode = regexec(&snapshotRegex, entry->d_name, 0, NULL, 0);
      if (retCode != 0) { continue; }

      std::string VmFilename = rootFilename + "/" + std::string(entry->d_name) + "/" + VmSubfile + "#";

      //Do we have a Vm file present?
      if (access((VmFilename + "000000").c_str(), R_OK) == -1) { continue; }

      std::shared_ptr<ParGridFunction> gf_Vm;
      double time;
      {
         //Open and read the file
         PFILE* file = Popen(VmFilename.c_str(), "r", COMM_LOCAL);
	 gf_Vm = std::make_shared<mfem::ParGridFunction>(pfespace);
    
	 OBJECT* hObj = file->headerObject;
	 std::vector<std::string> fieldNames,  fieldTypes;
	 objectGetv(hObj, "field_names", fieldNames); // field_names = gid Vm dVmD dVmR;
	 objectGetv(hObj, "field_types", fieldTypes); // field_types = u f f f;
         objectGet(hObj, "time", time, "-1");
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
      }
   

      StartTimer("Solve");
      
      heart_mat.Mult(*gf_Vm, gf_b);
      a->FormLinearSystem(ess_tdof_list,gf_x,gf_b,torso_mat,phi_e,phi_b);


      //HypreSmoother M(torso_mat, 6 /*GS*/);
      // PCG(torso_mat, M, phi_b, phi_e, 1, 2000, 1e-12);//, 0.0);
      
      pcg.Mult(phi_b,phi_e);

      EndTimer();
      
      // 11. Recover the solution as a finite element grid function.
      a->RecoverFEMSolution(phi_e, phi_b, gf_x);

      // 12. Save the refined mesh and the solution. This output can be viewed later
      //     using GLVis: "glvis -m refined.mesh -g sol.gf".
      std::ofstream sol_ofs(std::string(entry->d_name) + "/sol.gf");

      sol_ofs.precision(8);
      gf_x.Save(sol_ofs);

      for (int ielec=0; ielec<fileFromElectrode.size(); ielec++)
      {
         fileFromElectrode[ielec] << time << "\t" << gf_x[gfidFromElectrode[ielec]] << std::endl;
      }      
   }

   regfree(&snapshotRegex);
   closedir(dir);

   // Not sure how this will adapt to n>1, probably want to only run from rank 0 at least
   std::ofstream mesh_ofs("refined.mesh");
   mesh_ofs.precision(8);
   pmesh->Print(mesh_ofs);

   // 14. Free the used memory.
   delete M_test;
   delete a;
   delete b;
   delete pfespace;
   if (order > 0) { delete fec; }
   delete pmesh;
   delete mesh;

   return 0;
}
