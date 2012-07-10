#include "assignCellsToTasks.hh"
#include <string>
#include <cassert>
#include <cstdlib>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "object_cc.hh"
#include "ioUtils.h"
#include "Koradi.hh"
#include "writeCells.hh"
#include "Simulate.hh"
#include "GDLoadBalancer.hh"
#include "mpiUtils.h"
#include "GridPoint.hh"
#include "PerformanceTimers.hh"
#include "AnatomyCell.hh"

using namespace std;

namespace
{
    void koradiBalancer(Simulate& sim, OBJECT* obj, MPI_Comm comm);
    void gridBalancer(Simulate& sim, OBJECT* obj, MPI_Comm comm);
    void computeVolHistogram(Simulate& sim, vector<AnatomyCell>& cells, int nProcs, MPI_Comm comm);
}


void assignCellsToTasks(Simulate& sim, const string& name, MPI_Comm comm)
{
   OBJECT* obj = object_find(name.c_str(), "DECOMPOSITION");
   string method;
   objectGet(obj, "method", method, "koradi");

   profileStart("Assignment");
   if (method == "koradi")
      koradiBalancer(sim, obj, comm);
   else if (method == "grid")
      gridBalancer(sim, obj, comm);
   else
      assert(1==0);      
   profileStop("Assignment");
   sim.anatomy_.nLocal() = sim.anatomy_.cellArray().size();
   sim.anatomy_.nRemote() = 0;
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
      objectGet(obj, "outputRate",      kp.outputRate,      "-1");
      objectGet(obj, "alpha",           kp.alpha,           "0.05");
      objectGet(obj, "tolerance",       kp.tolerance,       "0.01");
      objectGet(obj, "nbrDeltaR",       kp.nbrDeltaR,       "2");
      
      kp.nCentersPerTask = nCenters/nTasks;
      assert(nCenters > 0);
      
      profileStart("Koradi");
      Koradi balancer(sim.anatomy_, kp);
      profileStop("Koradi");

      vector<AnatomyCell>& cells = sim.anatomy_.cellArray();
      computeVolHistogram(sim,cells,nCenters,comm);

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

      if (myRank == 0)
         cout << "GDLoadBalancer:  global grid " << sim.nx_ << " x " << sim.ny_ << " x " << sim.nz_ << ", process grid " << npex << " x " << npey << " x " << npez << endl;
      
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

      // allow user to define reduced process grid
      int rnx,rny,rnz;
      objectGet(obj, "rnx", rnx, "-1");
      objectGet(obj, "rny", rny, "-1");
      objectGet(obj, "rnz", rnz, "-1");
      if ((rnx != npex && rnx > 0) && (rny != npey && rny > 0) && (rnz != npez && rnz > 0))
         loadbal.setReducedProcGrid(rnx,rny,rnz);
      
      // compute initial data decomposition
      profileStart("gd_assign_init");
      loadbal.initialDistByVol(cells, sim.nx_, sim.ny_, sim.nz_);
      profileStop("gd_assign_init");

      computeVolHistogram(sim,cells,npegrid,comm);

      // set up visualization of initial distribution
      int visgrid;
      objectGet(obj, "visgrid", visgrid, "0");
      if (visgrid == 1)
      {
        stringstream name;
        name << "snap.gridinit";
        string fullname = name.str();
        if (myRank == 0)
           DirTestCreate(fullname.c_str());
        fullname += "/anatomy";
        writeCells(sim.anatomy_.cellArray(), sim.nx_, sim.ny_, sim.nz_, fullname.c_str());
      }
      
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

namespace
{
    void computeVolHistogram(Simulate& sim, vector<AnatomyCell>& cells, int nProcs, MPI_Comm comm)
    {
       int nTasks, myRank;
       MPI_Comm_size(comm, &nTasks);
       MPI_Comm_rank(comm, &myRank);
       bool testingOnly = (nProcs != nTasks);       
       
       // compute bounding box volumes from cells array
       vector<int> peminx(nProcs,99999999);
       vector<int> peminy(nProcs,99999999);
       vector<int> peminz(nProcs,99999999);
       vector<int> pemaxx(nProcs,-99999999);
       vector<int> pemaxy(nProcs,-99999999);
       vector<int> pemaxz(nProcs,-99999999);
       vector<int> nloc(nProcs,0);
       for (unsigned ii=0; ii<cells.size(); ++ii)
       {
          int peid = cells[ii].dest_;

          GridPoint gpt(cells[ii].gid_,sim.nx_,sim.ny_,sim.nz_);
          if (gpt.x < peminx[peid]) peminx[peid] = gpt.x;
          if (gpt.y < peminy[peid]) peminy[peid] = gpt.y;
          if (gpt.z < peminz[peid]) peminz[peid] = gpt.z;
          if (gpt.x > pemaxx[peid]) pemaxx[peid] = gpt.x;
          if (gpt.y > pemaxy[peid]) pemaxy[peid] = gpt.y;
          if (gpt.z > pemaxz[peid]) pemaxz[peid] = gpt.z;
          nloc[peid]++;
       }
       vector<int> peminx_all(nProcs);
       vector<int> peminy_all(nProcs);
       vector<int> peminz_all(nProcs);
       vector<int> pemaxx_all(nProcs);
       vector<int> pemaxy_all(nProcs);
       vector<int> pemaxz_all(nProcs);
       vector<int> nall(nProcs);
       MPI_Allreduce(&peminx[0], &peminx_all[0], nProcs, MPI_INT, MPI_MIN, comm);
       MPI_Allreduce(&peminy[0], &peminy_all[0], nProcs, MPI_INT, MPI_MIN, comm);
       MPI_Allreduce(&peminz[0], &peminz_all[0], nProcs, MPI_INT, MPI_MIN, comm);
       MPI_Allreduce(&pemaxx[0], &pemaxx_all[0], nProcs, MPI_INT, MPI_MAX, comm);
       MPI_Allreduce(&pemaxy[0], &pemaxy_all[0], nProcs, MPI_INT, MPI_MAX, comm);
       MPI_Allreduce(&pemaxz[0], &pemaxz_all[0], nProcs, MPI_INT, MPI_MAX, comm);
       MPI_Allreduce(&nloc[0], &nall[0], nProcs, MPI_INT, MPI_SUM, comm);

      if (myRank == 0)
      {
         const int nhistmax = 100; // number of bins
         vector<int> phist(nhistmax,0);
         int maxvol = 0;
         int minvol = sim.nx_*sim.ny_*sim.nz_;
         int maxvolip;
      
         vector<int> pevol_all(nProcs);
         for (int ip=0; ip<nProcs; ip++)
         {
            int tvol = 0;
            if (nall[ip] > 0)
               tvol = (pemaxx_all[ip]-peminx_all[ip]+1)*(pemaxy_all[ip]-peminy_all[ip]+1)*(pemaxz_all[ip]-peminz_all[ip]+1);
            if (tvol > maxvol)
            {
               maxvol = tvol;
               maxvolip = ip;
            }
            if (tvol < minvol) minvol = tvol;
            pevol_all[ip] = tvol;
         }
      
         int nhist = maxvol - minvol + 1;
         if (nhist > nhistmax) nhist = nhistmax;
         int delta = (maxvol-minvol + 1)/nhist;
         if ((maxvol-minvol+1)%nhist !=0) delta++;
         for (int ip=0; ip<nProcs; ip++)
         {
            int pvol = pevol_all[ip];
            int bin = (pvol-minvol)/delta;
            phist[bin]++;
         }
         cout << "load balance histogram (volume):  " << endl;
         for (int i=0; i<nhist; i++)
            cout << "  " << minvol+delta*i << " - " << minvol+delta*(i+1) << ":    " << phist[i] << endl;

         int voltot = 0;
         for (unsigned ip=0; ip<nProcs; ++ip)
            voltot += pevol_all[ip];
         double volavg = (double)voltot/(double)nProcs; 
         cout << "total assigned volume = " << voltot << ", avg. volume = " << volavg << ", max volume = " << maxvol << " (pe " << maxvolip << ")" << endl << endl;
      }
    }



}
