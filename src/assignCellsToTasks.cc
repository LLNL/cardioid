#include "assignCellsToTasks.hh"
#include <string>
#include <cassert>
#include <sstream>
#include <iostream>
#include <iomanip>

#include "object_cc.hh"
#include "ioUtils.h"
#include "Koradi.hh"
#include "writeCells.hh"
#include "Simulate.hh"
#include "GDLoadBalancer.hh"
#include "mpiUtils.h"
#include "GridPoint.hh"

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
      
      KoradiParms kp;
      int nCenters;      
      objectGet(obj, "nCenters", nCenters, buf.str());
      objectGet(obj, "verbose",         kp.verbose,         "1");
      objectGet(obj, "maxVoronoiSteps", kp.maxVoronoiSteps, "50");
      objectGet(obj, "maxSteps",        kp.maxSteps,        "500");
      objectGet(obj, "outputRate",      kp.outputRate,      "5");
      objectGet(obj, "alpha",           kp.alpha,           "0.05");
      objectGet(obj, "tolerance",       kp.tolerance,       "0.01");
      objectGet(obj, "nbrDeltaR",       kp.nbrDeltaR,       "2");
      
      kp.nCentersPerTask = nCenters/nTasks;
      assert(nCenters > 0);
      
      Koradi balancer(sim.anatomy_, kp);
      
      if (kp.nCentersPerTask > 1)
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

      // get process grid info 
      int npex, npey, npez;
      objectGet(obj, "nx", npex, "0");
      objectGet(obj, "ny", npey, "0");
      objectGet(obj, "nz", npez, "0");
      int npegrid = npex*npey*npez;
      assert(npegrid > 0);

      // Move each cell of gid=(i,j,k) to corresponding rank k
      for (unsigned ii=0; ii<cells.size(); ++ii)
      {
        GridPoint gpt(cells[ii].gid_,sim.nx_,sim.ny_,sim.nz_);
        if (nTasks >= sim.nz_) 
          cells[ii].dest_ = gpt.z;
        else
          cells[ii].dest_ = gpt.z*nTasks/sim.nz_;
      }

      sort(cells.begin(),cells.end(),AnatomyCell::destLessThan);
      unsigned nLocal = cells.size();
      vector<unsigned> dest(nLocal);
      for (unsigned ii=0; ii<cells.size(); ++ii)
         dest[ii] = cells[ii].dest_;

      // compute maximum possible number of non-zero local elements
      int nMax = sim.nx_*sim.ny_; 
      if (nTasks < sim.nz_)
        nMax *= sim.nz_/nTasks + 1;

      // carry out communication
      cells.resize(nMax);
      assignArray((unsigned char*)&(cells[0]), &nLocal, cells.capacity(),
                  sizeof(AnatomyCell), &(dest[0]), 0, comm);
      assert(nLocal <= nMax);
      cells.resize(nLocal);

      GDLoadBalancer loadbal(comm, npex, npey, npez);
      
      // compute initial data decomposition
      sim.tmap_["gd_assign_init"].start();
      loadbal.initialDistribution(cells, sim.nx_, sim.ny_, sim.nz_);
      sim.tmap_["gd_assign_init"].stop();

      //ewd DEBUG set up visualization of initial distribution
      int visgrid;
      objectGet(obj, "visgrid", visgrid, "0");
      if (visgrid == 1)
      {
        stringstream name;
        name << "snap.gridinit";
        string fullname = name.str();
        DirTestCreate(fullname.c_str());
        fullname += "/anatomy";
        writeCells(sim.anatomy_.cellArray(), sim.nx_, sim.ny_, sim.nz_, fullname.c_str());
      }
      
      // redistribute data until load balance is achieved
      int ninner;
      int threshold;
      objectGet(obj, "ninner", ninner, "10");
      objectGet(obj, "threshold", threshold, "4");
      int nmax = 10000000;

      sim.tmap_["gd_balance"].start();
      loadbal.balanceLoop(cells,ninner,threshold,nmax);
      sim.tmap_["gd_balance"].stop();

      //ewd DEBUG:  print out load balance data
      if (visgrid == 1)
      {
        stringstream name;
        name << "snap.gridfinal";
        string fullname = name.str();
        DirTestCreate(fullname.c_str());
        fullname += "/anatomy";
        writeCells(sim.anatomy_.cellArray(), sim.nx_, sim.ny_, sim.nz_, fullname.c_str());

        //ewd DEBUG
        //for (unsigned ii=0; ii<cells.size(); ++ii)
        //{
        //  if (cells[ii].dest_ == 100)
        //  {
        //    GridPoint gpt(cells[ii].gid_,sim.nx_,sim.ny_,sim.nz_);
        //    cout << " X  " << gpt.x << "  " << gpt.y << "  " << gpt.z << endl;
        //  }
        //}
        //ewd DEBUG

        
      }
      //ewd DEBUG
      
      //ewd print out load balance timing
      /*
      if (myRank == 0)
      {
        cout.setf(ios::fixed,ios::floatfield);
        for ( map<std::string,Timer>::iterator i = sim.tmap_.begin(); i != sim.tmap_.end(); i++ ) {
          double time = (*i).second.real();
          cout << "timing name=" << setw(15) << (*i).first << ":   time=" << setprecision(6) << setw(12) << time << " sec" << endl;
        }
      }
      */
      // ewd DEBUG


      bool testingOnly = (npegrid != nTasks);
      if (testingOnly)
      {
        if (myRank == 0)
          cout << "Process grid size does not match number of tasks.\n"
               << "This must be a test run only.  Exiting now"<< endl;
        exit(0);
      }
   }
}

