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
#include "BlockLoadBalancer.hh"
#include "workBoundBalancer.hh"
#include "pioBalancer.hh"
#include "mpiUtils.h"
#include "GridPoint.hh"
#include "PerformanceTimers.hh"
#include "AnatomyCell.hh"

using namespace std;

namespace
{
   void koradiBalancer(Simulate& sim, OBJECT* obj, MPI_Comm comm);
   void gridBalancer(Simulate& sim, OBJECT* obj, MPI_Comm comm);
   void blockBalancer(Simulate& sim, OBJECT* obj, MPI_Comm comm);
   int  workBoundScan(Simulate& sim, OBJECT* obj, MPI_Comm comm);
   int  pioBalancerScan(Simulate& sim, OBJECT* obj, MPI_Comm comm);

   void computeVolHistogram(Simulate& sim, vector<AnatomyCell>& cells, int nProcs, MPI_Comm comm);
   void computeNCellsHistogram(Simulate& sim, vector<AnatomyCell>& cells, int nProcs, MPI_Comm comm);
}


void assignCellsToTasks(Simulate& sim, const string& name, MPI_Comm comm)
{
   OBJECT* obj = object_find(name.c_str(), "DECOMPOSITION");
   string method;
   objectGet(obj, "method", method, "koradi");

   int nDiffusionCores=-1; 
   profileStart("Assignment");
   if (method == "koradi")
      koradiBalancer(sim, obj, comm);
   else if (method == "grid")
      gridBalancer(sim, obj, comm);
   else if (method == "block")
      blockBalancer(sim, obj, comm);
   else if (method == "workBound")
      nDiffusionCores = workBoundScan(sim, obj, comm);
   else if (method == "pio")
      nDiffusionCores = pioBalancerScan(sim, obj, comm);
   else
      assert(1==0);      
   profileStop("Assignment");
   sim.anatomy_.nLocal() = sim.anatomy_.cellArray().size();
   sim.anatomy_.nRemote() = 0;
}

namespace
{
   int  workBoundScan(Simulate& sim, OBJECT* obj, MPI_Comm comm)
   {
      int nTasks, myRank;
      MPI_Comm_size(comm, &nTasks);
      MPI_Comm_rank(comm, &myRank);

      vector<AnatomyCell>& cells = sim.anatomy_.cellArray();

      // get block dimension
      int dx,dy,dz;
      int target; 
      int printStats; 
      char defaultTarget[16];
      sprintf(defaultTarget,"%d",nTasks); 
      objectGet(obj, "dx", dx, "0");
      objectGet(obj, "dy", dy, "0");
      objectGet(obj, "dz", dz, "0");
      objectGet(obj, "printStats", printStats, "0");
      objectGet(obj, "targetNTasks",target,defaultTarget); 
      assert(dx*dy*dz != 0);
      assert((dz+2)%4 ==0); 
      int nx = sim.anatomy_.nx(); 
      int ny = sim.anatomy_.ny(); 
      int nz = sim.anatomy_.nz(); 
      int nDiffusionCores=workBoundBalancer(cells,dx,dy,dz,nx,ny,nz,target,printStats,comm); 

      computeNCellsHistogram(sim,cells,target,comm);
      computeVolHistogram(sim,cells,target,comm);

      if (target != nTasks) exit(0); 
      return nDiffusionCores; 
   }
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
      computeNCellsHistogram(sim,cells,nCenters,comm);
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

      computeNCellsHistogram(sim,cells,npegrid,comm);
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

    // load balancer uses Jim's algorithm which breaks work into uniform blocks
    // for dense regions, larger blocks for diffuse regions sized according to work
    // function
    
   void blockBalancer(Simulate& sim, OBJECT* obj, MPI_Comm comm)
   {
      int nTasks, myRank;
      MPI_Comm_size(comm, &nTasks);
      MPI_Comm_rank(comm, &myRank);

      vector<AnatomyCell>& cells = sim.anatomy_.cellArray();

      // get block dimension
      int bx,by,bz;
      objectGet(obj, "blockx", bx, "0");
      objectGet(obj, "blocky", by, "0");
      objectGet(obj, "blockz", bz, "0");
      assert(bx*by*bz != 0);
      
      // distribute columns to cells to processors:  column (i,j) will be on pe (i + ncx*j)
      int ncx = (sim.nx_%bx == 0 ? sim.nx_/bx : sim.nx_/bx + 1);
      int ncy = (sim.ny_%by == 0 ? sim.ny_/by : sim.ny_/by + 1);
      int nCols = ncx*ncy;
      
      int xmult = 1;
      int ymult = 1;
      bool toggle = false;
      while (nCols > xmult*ymult*nTasks)
      {
         if (toggle)
         {
            xmult *= 2;
            toggle = false;
         }
         else
         {
            ymult *= 2;
            toggle = true;
         }
      }
      
      int rcx = (ncx%xmult == 0 ? ncx/xmult : ncx/xmult + 1);
      int rcy = (ncy%ymult == 0 ? ncy/ymult : ncy/ymult + 1);
      int redTasks = rcx*rcy;
      
      if (myRank == 0 && xmult*ymult > 1)
         cout << "Block Load Balancer:  nCols = " << nCols << " (" << ncx << " x " << ncy << ")" << " distributed on " << redTasks << " tasks (" << rcx << " x " << rcy << ")" << endl;

      vector<int> dataOwner(nCols);      
      vector<int> localCols;
      for (unsigned icx=0; icx<ncx; ++icx)
      {
         for (unsigned icy=0; icy<ncy; ++icy)
         {
            int icol = icx + ncx*icy;
            int ipe = icx/xmult + rcx*icy/ymult;
            dataOwner[icol] = ipe;
            if (ipe == myRank)
               localCols.push_back(icol);
         }
      }

      // compute maximum possible size of cells array, carry out initial redistribution
      {
         vector<int> peCount(nTasks,0);
         for (unsigned ii=0; ii<cells.size(); ++ii)
         {
            GridPoint gpt(cells[ii].gid_,sim.nx_,sim.ny_,sim.nz_);
            int icol = (gpt.x/bx) + ncx*(gpt.y/by);
            int tdest = dataOwner[icol];
            cells[ii].dest_ = tdest;
            peCount[tdest]++;
         }
         vector<int> peCntSum(nTasks);
         MPI_Allreduce(&peCount[0], &peCntSum[0], nTasks, MPI_INT, MPI_SUM, comm);
         
         int nMax = 0;
         for (int ii=0; ii<nTasks; ii++)
            if (peCntSum[ii] > nMax) nMax = peCntSum[ii];
         
         sort(cells.begin(),cells.end(),AnatomyCell::destLessThan);
         unsigned nLocal = cells.size();
         vector<unsigned> dest(nLocal);
         for (unsigned ii=0; ii<cells.size(); ++ii)
            dest[ii] = cells[ii].dest_;
         
         // compute maximum possible number of non-zero local elements
         //      int nMax = xmult*ymult*bx*by*sim.nz_; 
         
         // carry out communication
         cells.resize(nMax);
         assignArray((unsigned char*)&(cells[0]), &nLocal, cells.capacity(),
                     sizeof(AnatomyCell), &(dest[0]), 0, comm);
         assert(nLocal <= nMax);
         cells.resize(nLocal);
      }

      BlockLoadBalancer loadbal(comm, sim.nx_, sim.ny_, sim.nz_, bx, by, bz);
      
      // compute decomposition using user-supplied parameters
      // (once we know this works, we'll replace with a tuning loop
      // to iterate parameters until load distributed across nTasks)
      double diffCost;
      objectGet(obj, "diffCost", diffCost, "0.125");
      if (myRank == 0)
         cout << "Block load balance starting with diffusion cost = " << diffCost << endl;

      profileStart("block_loadbal");
      int npes;
      int pesum_loc = 0;
      pesum_loc += loadbal.block(cells, diffCost, nCols, localCols);
      MPI_Allreduce(&pesum_loc, &npes, 1, MPI_INT, MPI_SUM, comm);
      if (myRank == 0)
            cout << "Diffusion cost = " << diffCost << ", load distributed on " << npes << " tasks." << endl;

      int optimizeTasks;
      objectGet(obj, "optimize", optimizeTasks, "0");
      int targetTasks;
      objectGet(obj, "nTasks", targetTasks, "-1");
      if (targetTasks == -1)
         targetTasks = nTasks;

      int biter = 0;
      const int biterMax = 100;
      bool optFound = false;
      double trialDiffCost = diffCost;
      bool tooMany = true;
      if (npes < targetTasks)
         tooMany = false;
      
      while (optimizeTasks && biter < biterMax && !optFound)
      {
         if (tooMany)
            trialDiffCost *= 0.99;
         else
            //trialDiffCost *= 1.01;
            trialDiffCost += 0.01;
            
         int pesum_loc = 0;
         pesum_loc += loadbal.block(cells, trialDiffCost, nCols, localCols);
         
         MPI_Allreduce(&pesum_loc, &npes, 1, MPI_INT, MPI_SUM, comm);
         if (myRank == 0)
            cout << "Diffusion cost = " << trialDiffCost << ", load distributed on " << npes << " tasks." << endl;
      }
      profileStop("block_loadbal");

      
      // collect task information from all tasks to work can be assigned to unique
      // MPI tasks
      vector<int> taskLoc(nTasks,0);
      vector<int> taskDist(nTasks,0);
      taskLoc[myRank] = pesum_loc;
      MPI_Allreduce(&taskLoc[0], &taskDist[0], nTasks, MPI_INT, MPI_SUM, comm);

      int peoff = 0;
      for (int ii=0; ii<myRank; ii++)
         peoff += taskDist[ii];

      for (unsigned ii=0; ii<cells.size(); ++ii)
      {
         int peloc = cells[ii].dest_;
         cells[ii].dest_ = peloc + peoff;
      }
      
      computeNCellsHistogram(sim,cells,npes,comm);
      computeVolHistogram(sim,cells,npes,comm);

      // output visualization data
      int printvis;
      objectGet(obj, "printvis", printvis, "0");
      if (printvis == 1)
      {
        stringstream name;
        name << "block.loadbal";
        string fullname = name.str();
        if (myRank == 0)
           DirTestCreate(fullname.c_str());
        fullname += "/anatomy";
        writeCells(sim.anatomy_.cellArray(), sim.nx_, sim.ny_, sim.nz_, fullname.c_str());
      }
      
      bool testingOnly = (npes > nTasks);
      if (testingOnly)
      {
         if (myRank == 0)
         {
            cout << "Load distributed over " << npes << " tasks, greater than available " << nTasks << " tasks." << endl;
            cout << "This must be a test run only.  Exiting now"<< endl;
         }
         exit(0);
      }
      else
      {
         // compute maximum possible size of cells array, carry out final redistribution
         vector<int> peCount(nTasks,0);
         for (unsigned ii=0; ii<cells.size(); ++ii)
         {
            int tdest = cells[ii].dest_;
            peCount[tdest]++;
         }
         vector<int> peCntSum(nTasks);
         MPI_Allreduce(&peCount[0], &peCntSum[0], nTasks, MPI_INT, MPI_SUM, comm);
         
         int nMax = 0;
         for (int ii=0; ii<nTasks; ii++)
            if (peCntSum[ii] > nMax) nMax = peCntSum[ii];
         
         sort(cells.begin(),cells.end(),AnatomyCell::destLessThan);
         unsigned nLocal = cells.size();
         vector<unsigned> dest(nLocal);
         for (unsigned ii=0; ii<cells.size(); ++ii)
            dest[ii] = cells[ii].dest_;
         
         // compute maximum possible number of non-zero local elements
         //      int nMax = xmult*ymult*bx*by*sim.nz_; 
         
         // carry out communication
         cells.resize(nMax);
         assignArray((unsigned char*)&(cells[0]), &nLocal, cells.capacity(),
                     sizeof(AnatomyCell), &(dest[0]), 0, comm);
         assert(nLocal <= nMax);
         cells.resize(nLocal);
      }
   }
}


namespace
{
   int pioBalancerScan(Simulate& sim, OBJECT* obj, MPI_Comm comm)
   {
      string domainFile;
      string pxyzFile;
      objectGet(obj, "domainFile", domainFile, "balance/domains#");
      objectGet(obj, "pxyzFile",   pxyzFile,   "balance/pxyz#");

      int nD = pioBalancer(domainFile, pxyzFile, sim, comm);
      return nD;
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

namespace
{
    void computeNCellsHistogram(Simulate& sim, vector<AnatomyCell>& cells, int nProcs, MPI_Comm comm)
    {
       int nTasks, myRank;
       MPI_Comm_size(comm, &nTasks);
       MPI_Comm_rank(comm, &myRank);
       bool testingOnly = (nProcs != nTasks);       
       
       // compute load histogram from data in cells
       vector<int> histnloc_(nProcs,0);
       vector<int> mydata(nProcs,0);
       for (unsigned ii=0; ii<cells.size(); ++ii)
          mydata[cells[ii].dest_]++;
       MPI_Allreduce(&mydata[0], &histnloc_[0], nProcs, MPI_INT, MPI_SUM, comm);
       if (myRank == 0)
       {
          // compute histogram of current data distribution
          const int nhistmax = 100; // number of bins
          vector<int> phist(nhistmax,0);

          int maxnum = 0;
          int minnum = sim.nx_*sim.ny_*sim.nz_;
          int maxpe = -1;
          for (int p=0; p<nProcs; p++)
          {
             if (histnloc_[p] > maxnum) {
                maxnum = histnloc_[p];
                maxpe = p;
             }
             if (histnloc_[p] < minnum)
                minnum = histnloc_[p];
          }
          int nhist = maxnum - minnum + 1;
          if (nhist > nhistmax) nhist = nhistmax;
          int delta = (maxnum-minnum + 1)/nhist;
          if ((maxnum-minnum+1)%nhist !=0) delta++;
          for (int p=0; p<nProcs; p++)
          {
             int bin = (histnloc_[p]-minnum)/delta;
             phist[bin]++;
          }
          cout << "load balance histogram (ncells):  " << endl;
          for (int i=0; i<nhist; i++)
             cout << "  " << minnum+delta*i << " - " << minnum+delta*(i+1) << ":    " << phist[i] << endl;
          
          int nloctot_ = 0;
          for (unsigned ii=0; ii<nProcs; ++ii)
             nloctot_ += histnloc_[ii];
          double nlocavg_ = (double)nloctot_/(double)nProcs; 
          
          cout << "total # of non-zero grid points = " << nloctot_ << ", avg. # per task = " << nlocavg_ << ", max pe = " << maxpe << " (" << maxnum << " cells)" << endl << endl;

       }
    }
}
