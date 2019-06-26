#include "mfem.hpp"
#include "object.h"
#include "object_cc.hh"
#include "ddcMalloc.h"
#include "pio.h"
#include "pioFixedRecordHelper.h"
#include "units.h"
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
#include <sys/stat.h>
#include "util.hpp"
#include "MatrixElementPiecewiseCoefficient.hpp"
#include "cardiac_coefficients.hpp"

#define StartTimer(x)
#define EndTimer()

using namespace mfem;

MPI_Comm COMM_LOCAL = MPI_COMM_WORLD;

//Stolen from SingleCell
class Timeline
{
 public:
   Timeline(double dt, double duration)
   {
      maxTimesteps_ = round(duration/dt);
      dt_ = duration/maxTimesteps_;
   }
   int maxTimesteps() const { return maxTimesteps_; };
   double dt() const { return dt_; }
   double maxTime() const { return dt_*maxTimesteps_; }
   double realTimeFromTimestep(int timestep) const
   {
      return timestep*dt_;
   }
   int timestepFromRealTime(double realTime) const
   {
      return round(realTime/dt_);
   }
   std::string outputIdFromTimestep(const int timestep) const
   {
      double resolution = 1e-3;
      int width = 8;
      while (resolution > dt_) {
         resolution /= 10;
         width++;
      }
      std::stringstream ss;
      ss << std::setfill('0') << std::setw(width)
         << int(round(dt_*timestep/resolution));
      return ss.str();
   }

 private:
   double dt_;
   int maxTimesteps_;
};

class OutputCoordinator
{

 private:
   
};


void recursive_mkdir(const std::string dirname, mode_t mode=S_IRWXU|S_IRWXG)
{
   int startSearch=0;
   do
   {
      int endSearch = dirname.find("/", startSearch);
      //if directory doesn't exist
      if (endSearch < 0) {
         endSearch = dirname.length();
      }
      std::string thisDirname = dirname.substr(0, endSearch);
      DIR* dir = opendir(thisDirname.c_str());
      if (dir)
      {
         closedir(dir);
      }
      else if (ENOENT == errno) {
         //make the directory
         int ret = mkdir(thisDirname.c_str(), mode);
         assert(ret == 0);
      }
      startSearch=endSearch+1;
   } while (startSearch < dirname.length());
}


void save1dNumpyArray(const std::string filename, const std::vector<double>& data)
{
   const int NUMPY_HEADER_SIZE = 128;
   char numpyHeaderFull[NUMPY_HEADER_SIZE];
   //FIXME, report the correct endianness.  I'm too lazy ATM.
   char numpyHeader1[] = 
      "\x93NUMPY\x01\x00\x76\x00{'descr': '<f8', 'fortran_order': False, 'shape': (";
   char numpyHeader2[] = ",)}";
   int cursor=0;
   memcpy(numpyHeaderFull+cursor, numpyHeader1, sizeof(numpyHeader1)-1);
   cursor += sizeof(numpyHeader1)-1;
   {
      std::stringstream ss;
      ss << data.size();
      ss.str();
      memcpy(numpyHeaderFull+cursor,ss.str().c_str(),ss.str().size());
      cursor += ss.str().size();
   }
   memcpy(numpyHeaderFull+cursor,numpyHeader2, sizeof(numpyHeader2)-1);
   cursor += sizeof(numpyHeader2)-1;
   for (; cursor<NUMPY_HEADER_SIZE-1; cursor++) {
      numpyHeaderFull[cursor] = ' ';
   }
   numpyHeaderFull[NUMPY_HEADER_SIZE-1] = '\n';
   FILE* numpyFile = fopen(filename.c_str(), "w");
   assert(numpyFile != NULL);
   std::size_t bytesWritten = fwrite(numpyHeaderFull, sizeof(char), NUMPY_HEADER_SIZE, numpyFile);
   assert(bytesWritten == NUMPY_HEADER_SIZE);
   std::size_t doublesWritten = fwrite(&data[0], sizeof(double), data.size(), numpyFile);
   assert(doublesWritten == data.size());
   fclose(numpyFile);
}

int main(int argc, char *argv[])
{
   MPI_Init(NULL,NULL);
   int num_ranks, my_rank;
   MPI_Comm_size(COMM_LOCAL,&num_ranks);
   MPI_Comm_rank(COMM_LOCAL,&my_rank);

   units_internal(1e-3, 1e-9, 1e-3, 1e-3, 1, 1e-9, 1);
   units_external(1e-3, 1e-9, 1e-3, 1e-3, 1, 1e-9, 1);

   if (my_rank == 0)
   {
      std::cout << "Initializing with " << num_ranks << " MPI ranks." << std::endl;
   }
   
   int order = 1;

   std::vector<std::string> objectFilenames;
   if (argc == 1)
      objectFilenames.push_back("femheart.data");

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

   std::vector<double> sigma_m;
   objectGet(obj,"sigma_m",sigma_m);
   assert(heartRegions.size()*3 == sigma_m.size());

   double dt;
   objectGet(obj,"dt",dt,"0.01 ms");
   double Bm;
   objectGet(obj,"Bm",Bm,"140"); // 1/mm
   double Cm;
   objectGet(obj,"Cm",Cm,"0.01"); // 1 uF/cm^2 = 0.01 uF/mm^2
 
   std::string reactionName;
   objectGet(obj, "reaction", reactionName, "BetterTT06");

   std::string outputDir;
   objectGet(obj, "outdir", outputDir, ".");
   
   double endTime;
   objectGet(obj, "end_time", endTime, "0 ms");

   double outputRate;
   objectGet(obj, "output_rate", outputRate, "1 ms");

   //double checkpointRate;
   //objectGet(obj, "checkpoint_rate", checkpointRate, "100 ms");

   double initVm;
   objectGet(obj, "init_vm", initVm, "-83");

   bool useNodalIion;
   objectGet(obj, "nodal_ion", useNodalIion, "1");

   StimulusCollection stims(dt);
   {
      std::vector<std::string> stimulusNames;
      objectGet(obj, "stimulus", stimulusNames);
      for (auto name : stimulusNames)
      {
         OBJECT* stimobj = object_find(name.c_str(), "STIMULUS");
         assert(stimobj != NULL);
         int numTimes;
         objectGet(stimobj, "n", numTimes, "1");
         double bcl;
         objectGet(stimobj, "bcl", bcl, "0 ms");
         assert(numTimes == 1 || bcl != 0);
         double startTime;
         objectGet(stimobj, "start", startTime, "0 ms");
         double duration;
         objectGet(stimobj, "duration", duration, "1 ms");
         double strength;
         objectGet(stimobj, "strength", strength, "0"); //uA/uF
         assert(strength >= 0);
         std::string location;
         objectGet(stimobj, "where", location, "");
         assert(!location.empty());
         OBJECT* locobj = object_find(location.c_str(), "REGION");
         assert(locobj != NULL);
         std::string regionType;
         objectGet(locobj, "type", regionType, "");
         assert(!regionType.empty());
         shared_ptr<StimulusLocation> stimLoc;
         if (regionType == "ball")
         {
            std::vector<double> center;
            objectGet(locobj, "center", center);
            assert(center.size() == 3);
            double radius;
            objectGet(locobj, "radius", radius, "-1");
            assert(radius >= 0);
            stimLoc = std::make_shared<CenterBallStimulus>(center[0],center[1],center[2],radius);
         }
         else if (regionType == "box")
         {
            std::vector<double> lower;
            objectGet(locobj, "lower", lower);
            assert(lower.size() == 3);
            vector<double> upper;
            objectGet(locobj, "upper", upper);
            assert(upper.size() == 3);
            stimLoc = std::make_shared<BoxStimulus>
               (lower[0], upper[0],
                lower[1], upper[1],
                lower[2], upper[2]);
         }
         shared_ptr<StimulusWaveform> stimWave(new SquareWaveform());
         stims.add(Stimulus(numTimes, startTime, duration, bcl, strength, stimLoc, stimWave));
      }
   }
   
   Timeline timeline(dt, endTime);   

   /*
     Ok, we're solving:

     div(sigma_m*grad(Vm)) = Bm*Cm*(dVm/dt + Iion - Istim)

     time in ms
     space in mm
     Vm in mV
     Iion in uA/uF
     Istim in uA/uF
     Cm in uF/mm^2
     Bm in 1/mm
     sigma in mS/mm

     To solve this, I use a semi implicit crank nicolson:

     div(sigma_m*grad[(Vm_new+Vm_old)/2]) = Bm*Cm*[(Vm_new-Vm_old)/dt + Iion - Istim]

     with some algebra, that comes to

     {1+dt/(2*Bm*Cm)*(-div sigma_m*grad)}*Vm_new = {1-dt/(2*Bm*Cm)*(-div sigma_m*grad)}*Vm_old - dt*Iion + dt*Istim

     One easy way to check this is to set sigma_m to zero, then you get Forward euler for isolated equations.

     Each {} is a matrix that is done with FEM.

     sigma_m = sigma_e_diagonal*sigma_i_diagonal/(sigma_e_diagonal+sigma_i_diagonal)

     This is the monodomain conductivity.  It really only approximates the bidomain conductivity if the ratio
     of sigma_e_tensor = k*sigma_i_tensor.  When this happens, you can remove Phi_e from the equations and
     end up with the equation listed above. 
   */

   
   StartTimer("Setting Attributes");
   mesh->SetAttributes();
   EndTimer();

   StartTimer("Partition Mesh");
   // If I read correctly, pmeshpart will now point to an integer array
   //  containing a partition ID (rank!) for every element ID.
   int *pmeshpart = mesh->GeneratePartitioning(num_ranks);
   EndTimer();


   //Go through all the elements and label the partitioning for the vertices
   std::vector<set<int> > pvertset(mesh->GetNV());
   for (int ielem=0; ielem<mesh->GetNE(); ielem++)
   {
      Array<int> verts;
      mesh->GetElementVertices(ielem, verts);
      for (int ivert=0; ivert<verts.Size(); ivert++)
      {
         pvertset[verts[ivert]].insert(pmeshpart[ielem]);
      }
   }

   std::vector<int> local_extents(num_ranks+1);
   {
      std::vector<int> local_counts(num_ranks, 0);
      for(int i=0; i<mesh->GetNV(); i++)
      {
         if ( ! pvertset[i].empty())
         {
            local_counts[*(pvertset[i].begin())]++;
         }
      }

      local_extents[0] = 0;
      for (int irank=0; irank<num_ranks; irank++)
      {
         local_extents[irank+1] = local_extents[irank]+local_counts[irank];
      }
   }

   std::vector<int> globalvert_from_ranklookup(local_extents[num_ranks]);
   std::vector<int> ghostlocalvert_from_ranklookup(local_extents[num_ranks]);
   {
      std::vector<int> cursor_ghostlocal_from_rank(num_ranks, 0);
      std::vector<int> cursor_ranklookup_from_rank = local_extents;
      for(int i=0; i<mesh->GetNV(); i++)
      {
         if ( ! pvertset[i].empty())
         {
            int irank = *(pvertset[i].begin());
            int ranklookup = cursor_ranklookup_from_rank[irank]++;
            int globalvert = i;
            int ghostlocal = cursor_ghostlocal_from_rank[irank];
            globalvert_from_ranklookup[ranklookup] = globalvert;
            ghostlocalvert_from_ranklookup[ranklookup] = ghostlocal;
            for (const int used_by_this_rank : pvertset[i])
            {
               cursor_ghostlocal_from_rank[used_by_this_rank]++;
            }
         }
      }
   }

   //Get the element material types for each index.
   std::vector<int> material_from_ranklookup(local_extents[num_ranks]);
   {
      std::vector<int> element_from_globalvert(mesh->GetNV(), mesh->GetNE());
      for (int ielem=0; ielem<mesh->GetNE(); ielem++)
      {
         Array<int> verts;
         mesh->GetElementVertices(ielem, verts);
         for (int ivert=0; ivert<verts.Size(); ivert++)
         {
            element_from_globalvert[verts[ivert]] = std::min(element_from_globalvert[verts[ivert]], ielem);
         }
      }
      std::vector<int> cursor_ranklookup_from_rank = local_extents;
      for(int i=0; i<mesh->GetNV(); i++)
      {
         if ( ! pvertset[i].empty())
         {
            int irank = *(pvertset[i].begin());
            int ranklookup = cursor_ranklookup_from_rank[irank]++;
            int globalvert = i;

            int ielem = element_from_globalvert[globalvert];
            material_from_ranklookup[ranklookup] = mesh->GetElement(ielem)->GetAttribute();
         }
      }
   }
   
   if (my_rank == 0)
   {
      for(int i=0; i<num_ranks; i++) {
         std::cout << "Rank " << i << " has " << local_extents[i+1]-local_extents[i] << " nodes!" << std::endl;
      }
   }
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh, pmeshpart);
   
   // Build a new FEC...
   FiniteElementCollection *fec;
   if (my_rank == 0) { std::cout << "Creating new FEC..." << std::endl; }
   fec = new H1_FECollection(order, dim);
   // ...and corresponding FES
   ParFiniteElementSpace *pfespace = new ParFiniteElementSpace(pmesh, fec);
   FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);
   std::cout << "[" << my_rank << "] Number of finite element unknowns: "
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
   ParGridFunction gf_Vm(pfespace);
   ParGridFunction gf_b(pfespace);
   gf_Vm = initVm;
   gf_b = 0.0;

   // Load fiber quaternions from file
   std::shared_ptr<GridFunction> flat_fiber_quat;
   ecg_readGF(obj, "fibers", mesh, flat_fiber_quat);
   std::shared_ptr<ParGridFunction> fiber_quat;
   fiber_quat = std::make_shared<mfem::ParGridFunction>(pmesh, flat_fiber_quat.get(), pmeshpart);

   
   // Load conductivity data
   MatrixElementPiecewiseCoefficient sigma_m_pos_coeffs(fiber_quat);
   MatrixElementPiecewiseCoefficient sigma_m_neg_coeffs(fiber_quat);
   for (int ii=0; ii<heartRegions.size(); ii++) {
      int heartCursor=3*ii;
      Vector sigma_m_vec(&sigma_m[heartCursor],3);
      Vector sigma_m_pos_vec(3);
      Vector sigma_m_neg_vec(3);
      for (int jj=0; jj<3; jj++)
      {
         double value = sigma_m[jj]*dt/2/Bm/Cm;
         sigma_m_pos_vec[jj] = value;
         sigma_m_neg_vec[jj] = -value;
      }
    
      sigma_m_pos_coeffs.heartConductivities_[heartRegions[ii]] = sigma_m_pos_vec;
      sigma_m_neg_coeffs.heartConductivities_[heartRegions[ii]] = sigma_m_neg_vec;
   }

   // 8. Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the Laplacian operator -Delta, by adding the Diffusion
   //    domain integrator.

   // NOTICE THE FLIP IN SIGNS FOR SIGMA!  This is on purpose, Diffusion does -div(sigma*grad)

   StartTimer("Forming bilinear system (RHS)");

   ConstantCoefficient one(1.0);
   ParBilinearForm *b = new ParBilinearForm(pfespace);
   b->AddDomainIntegrator(new DiffusionIntegrator(sigma_m_neg_coeffs));
   b->AddDomainIntegrator(new MassIntegrator(one));
   b->Assemble();
   // This creates the linear algebra problem.
   HypreParMatrix RHS_mat;
   b->FormSystemMatrix(ess_tdof_list, RHS_mat);
   EndTimer();

   StartTimer("Forming bilinear system (LHS)");
   
   // Brought out of loop to avoid unnecessary duplication
   ParBilinearForm *a = new ParBilinearForm(pfespace);   // defines a.
   a->AddDomainIntegrator(new DiffusionIntegrator(sigma_m_pos_coeffs));
   a->AddDomainIntegrator(new MassIntegrator(one));
   a->Update(pfespace);
   a->Assemble();
   HypreParMatrix LHS_mat;
   a->FormSystemMatrix(ess_tdof_list,LHS_mat);
   EndTimer();

   //Set up the solve
   HyprePCG pcg(LHS_mat);
   pcg.SetTol(1e-12);
   pcg.SetMaxIter(2000);
   pcg.SetPrintLevel(2);
   HypreSolver *M_test = new HypreBoomerAMG(LHS_mat);
   pcg.SetPreconditioner(*M_test);


   //Set up the ionic models
   ParLinearForm *c = new ParLinearForm(pfespace);
   //positive dt here because the reaction models use dVm = -Iion
   c->AddDomainIntegrator(new DomainLFIntegrator(stims));


   
   ThreadServer& threadServer = ThreadServer::getInstance();
   ThreadTeam defaultGroup = threadServer.getThreadTeam(vector<unsigned>());
   std::vector<std::string> reactionNames;
   reactionNames.push_back(reactionName);
   std::vector<int> cellTypes;

   //int Iion_order = 2*order+3;
   int Iion_order = 2*order-1;
   QuadratureSpace quadSpace(pmesh, Iion_order);
   if (useNodalIion)
   {
      for (int ranklookup=local_extents[my_rank]; ranklookup<local_extents[my_rank+1]; ranklookup++)
      {
         cellTypes.push_back(material_from_ranklookup[ranklookup]);
      }
   }
   else
   {
      for (int i = 0; i < pfespace->GetNE(); ++i)
      {
         ElementTransformation *T = pfespace->GetElementTransformation(i);
         //This is a hack.  There's no way to get access to the offsets() array
         //in Quadrature Space without declaring ourselves to be a friend class.
         //This is broken and I hope it is fixed in 4.0
         Vector localVm;
         int NNN = quadSpace.GetElementIntRule(i).GetNPoints();
         for ( int j=0; j<NNN; j++)
         {
            cellTypes.push_back(T->Attribute);
         }
      }
   }

   ReactionWrapper reactionWrapper(dt,reactionNames,defaultGroup,cellTypes);
   reactionWrapper.Initialize();
   cellTypes.clear();
   reactionNames.clear();
   
   ParBilinearForm *Iion_blf;
   HypreParMatrix Iion_mat;
   ConstantCoefficient dt_coeff(dt);
   ReactionFunction* rf = NULL;
   if (useNodalIion) {
      Iion_blf = new ParBilinearForm(pfespace);
      Iion_blf->AddDomainIntegrator(new MassIntegrator(dt_coeff));
      Iion_blf->Update(pfespace);
      Iion_blf->Assemble();
      Iion_blf->FormSystemMatrix(ess_tdof_list,Iion_mat);
   } else {
      Iion_blf = NULL;
      
      rf = new ReactionFunction(&quadSpace,pfespace,&reactionWrapper); 
      c->AddDomainIntegrator(new QuadratureIntegrator(rf, dt)); 
   }

   Vector actual_Vm(pfespace->GetTrueVSize()), actual_b(pfespace->GetTrueVSize()), actual_old(pfespace->GetTrueVSize());
   Vector actual_Iion(pfespace->GetTrueVSize());
   bool first=true;

   if (useNodalIion)
   {
      actual_Vm = reactionWrapper.getVmReadonly();
   }
   
   int itime=0;
   while (1)
   {
      if (my_rank == 0)
      {
         std::cout << "time = " << timeline.realTimeFromTimestep(itime) << std::endl;
      }
      //output if appropriate
      if ((itime % timeline.timestepFromRealTime(outputRate)) == 0)
      {
         if (my_rank ==0)
         {
            std::string timedir = outputDir + "/tm" + timeline.outputIdFromTimestep(itime);
            recursive_mkdir(timedir); 

            std::vector<double> dataBuffer(local_extents[num_ranks]);
            for (int irank=0; irank<num_ranks; irank++)
            {
               int local_size = local_extents[irank+1] - local_extents[irank];
               std::vector<double> rankBuffer(local_size);
               if (irank==0)
               {
                  if (local_extents[irank+1] > 0)
                  {
                     //Since rank 0 is always the least, ranklookup == localvert
                     //the following assertion makes sure this is always true.
                     assert(ghostlocalvert_from_ranklookup[local_extents[irank+1]-1] == local_extents[irank+1]-1);
                     memcpy(&rankBuffer[0], &gf_Vm[0], sizeof(double)*local_size);
                  }
               }
               else
               {
                  MPI_Status dontcare;
                  MPI_Recv(&rankBuffer[0], local_size,
                           MPI_DOUBLE, irank, 455, MPI_COMM_WORLD, &dontcare);
               }
               for (int ii=0; ii<local_size; ii++)
               {
                  dataBuffer[globalvert_from_ranklookup[ii+local_extents[irank]]] = rankBuffer[ii];
               }
            }
            std::string VmFilename = timedir + "/Vm.npy";
            save1dNumpyArray(timedir + "/Vm.npy", dataBuffer);
         }
         else
         {
            int local_size = local_extents[my_rank+1]-local_extents[my_rank];
            std::vector<double> dataFromLocalRanklookup(local_size);
            for (int ii=0; ii<local_size; ii++)
            {
               int ranklookup = local_extents[my_rank] + ii;
               dataFromLocalRanklookup[ii] = gf_Vm[ghostlocalvert_from_ranklookup[ranklookup]];
            }
            MPI_Send(&dataFromLocalRanklookup[0], local_size,
                     MPI_DOUBLE, 0, 455, MPI_COMM_WORLD);
         }
      }
      //if end time, then exit
      if (itime == timeline.maxTimesteps()) { break; }

      //calculate the ionic contribution.
      if (useNodalIion) {
         reactionWrapper.getVmReadwrite() = actual_Vm; //should be a memcpy
         reactionWrapper.Calc();
      } else {
         rf->Calc(gf_Vm);
      }
      
      //add stimulii
      stims.updateTime(timeline.realTimeFromTimestep(itime));
      
      //compute the Iion and stimulus contribution
      c->Update();
      c->Assemble();
      a->FormLinearSystem(ess_tdof_list, gf_Vm, *c, LHS_mat, actual_Vm, actual_b, 1);
      //compute the RHS matrix contribution
      RHS_mat.Mult(actual_Vm, actual_old);
      actual_b += actual_old;

      if (useNodalIion)
      {
         Iion_mat.Mult(reactionWrapper.getIionReadonly(), actual_old);
         actual_b += actual_old;
      }
      //solve the matrix
      pcg.Mult(actual_b, actual_Vm);

      a->RecoverFEMSolution(actual_Vm, *c, gf_Vm);

      itime++;
      first=false;
   }

   // 14. Free the used memory.
   delete M_test;
   delete a;
   delete b;
   delete c;
   delete pfespace;
   if (order > 0) { delete fec; }
   delete mesh, pmesh, pmeshpart;
   
   return 0;
}
