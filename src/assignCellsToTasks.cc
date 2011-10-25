#include "assignCellsToTasks.hh"
#include <string>
#include <cassert>
#include <sstream>
#include <iostream>

#include "object_cc.hh"
#include "ioUtils.h"
#include "KoradiTest.hh"
#include "writeCells.hh"
#include "Simulate.hh"
#include "GDLoadBalancer.hh"
#include "mpiUtils.h"

using namespace std;

namespace
{
   void koradiBalancer(Simulate& sim, OBJECT* obj, MPI_Comm comm);
   void gridBalancer(Simulate& sim, OBJECT* obj, MPI_Comm comm);
}


void assignCellsToTasks(Simulate& sim, const string& name, MPI_Comm comm)
{
   OBJECT* obj = object_find(name.c_str(), "DECOMPOSITION");
   string method;
   objectGet(obj, "method", method, "koradi");

   if (method == "koradi")
      koradiBalancer(sim, obj, comm);
   else if (method == "grid")
      gridBalancer(sim, obj, comm);
   else
      assert(1==0);	 
}

namespace
{
   void koradiBalancer(Simulate& sim, OBJECT* obj, MPI_Comm comm)
   {
      int nTasks;  MPI_Comm_size(comm, &nTasks);
      stringstream buf;
      buf << nTasks;

      int nCenters;
      double alpha;
      int koradiSteps;
      int voronoiSteps;
      objectGet(obj, "nCenters", nCenters, buf.str());
      objectGet(obj, "koradiSteps", koradiSteps, "100");
      objectGet(obj, "voronoiSteps", voronoiSteps, "30");
      objectGet(obj, "alpha", alpha, "0.05");

      int nCentersPerTask = nCenters/nTasks;
      
      KoradiTest tester(sim, nCentersPerTask, alpha, voronoiSteps);
      for (unsigned ii=0; ii<koradiSteps; ++ii)
      {
	 tester.balanceStep();
	 stringstream name;
	 name << "snap."<<ii;
	 string fullname = name.str();
	 DirTestCreate(fullname.c_str());
	 fullname += "/anatomy";
	 writeCells(sim, fullname.c_str());

	 if (ii%100 == 0)
	    tester.recondition(voronoiSteps);

      }
      
      exit(0);
   }
   
}





namespace
{
   
   // GDLoadBalancer expects a linear array of cell types with size
   // nx*ny*nz.  For now, all cells are expected to be on rank 0.  The
   // types array is built from the data in Simulate, we call the
   // balancer, then use the resulting data to reassign the cells in the
   // Simulate anatomy.
   //
   // I think it would be not too hard to get the GDLoadBalancer to work
   // directly from the cells array, but I'll save that effort for Erik
   // (along with parallelizing the GDLoadBalancer).

   void gridBalancer(Simulate& sim, OBJECT* obj, MPI_Comm comm)
   {
      int nTasks, myRank;
      MPI_Comm_size(comm, &nTasks);
      MPI_Comm_rank(comm, &myRank);

      vector<AnatomyCell>& cells = sim.anatomy_.cellArray();
      
      // Move all cells to rank 0
      for (unsigned ii=0; ii<cells.size(); ++ii)
	 cells[ii].dest_ = 0;
      unsigned nLocal = cells.size();
      unsigned nGlobal;
      MPI_Allreduce(&nLocal, &nGlobal, 1, MPI_UNSIGNED, MPI_SUM, comm);
      vector<unsigned> dest(nLocal, 0);
      
      cells.resize(nGlobal);
      assignArray((unsigned char*)&(cells[0]), &nLocal, cells.capacity(),
		  sizeof(AnatomyCell), &(dest[0]), 0, comm);
      cells.resize(nLocal);
      
      // build the types array
      int nGrid = sim.nx_ * sim.ny_ * sim.nz_;
      vector<int> types(nGrid, 0);
      for (unsigned ii=0; ii<cells.size(); ++ii)
	 types[cells[ii].gid_] = cells[ii].cellType_;
      

      // get process grid info 
      int nx, ny, nz;
      objectGet(obj, "nx", nx, "0");
      objectGet(obj, "ny", ny, "0");
      objectGet(obj, "nz", nz, "0");
      assert(nx*ny*nz > 0);

      GDLoadBalancer loadbal(nx, ny, nz);
      
      // compute data decomposition, redistribute data until load
      // balance is achieved
      sim.tmap_["assign_init"].start();
      loadbal.initialDistribution(types, sim.nx_, sim.ny_, sim.nz_);
      sim.tmap_["assign_init"].stop();
      
      sim.tmap_["balance"].start();
      loadbal.balanceLoop();
      sim.tmap_["balance"].stop();

      // If we have more grid points than MPI ranks, this must be a test.
      if (nx*ny*nz != nTasks)
      {
	 if (myRank == 0)
	    cout << "Grid size does not match number of tasks.\n"
		 << "This must be a test run only.  Exiting now"<< endl;
	 exit(0);
      }

      // This is where we send cells to their tasks.
      // To make this in any way efficient we need to sort the cell
      // array so that we can match up with GDLoadBalancer::gpe_
      // I don't have the patience for that right now, and the whole
      // system would be worlds easier if GDLoadBalancer would just move
      // the cells around.  I'll negotiate with Erik later.

   }
}

