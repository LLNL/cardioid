#include "getRemoteCells.hh"
#include <vector>

#include "Simulate.hh"
#include "Anatomy.hh"
#include "GridRouter.hh"
#include "HaloExchange.hh"
#include "pio.h"
#include "object_cc.hh"

#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

void dumpCells(const Anatomy& anatomy)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, & myRank);
   {
      stringstream buf;
      buf << "cellsLocal."<<myRank;
      ofstream file(buf.str().c_str());
      for (unsigned ii=0; ii<anatomy.nLocal(); ++ii)
      {
         Tuple gg = anatomy.globalTuple(ii);
         file << gg.x() << " " << gg.y() << " " << gg.z() <<endl;
      }
   }
   {
      stringstream buf;
      buf << "cellsRemote."<<myRank;
      ofstream file(buf.str().c_str());
      for (unsigned ii=anatomy.nLocal(); ii<anatomy.size(); ++ii)
      {
         Tuple gg = anatomy.globalTuple(ii);
         file << gg.x() << " " << gg.y() << " " << gg.z() <<endl;
      }
   }
   
}



void getRemoteCells(Simulate& sim, const string& name, MPI_Comm comm)
{
   Anatomy& anatomy = sim.anatomy_;

   int nx = sim.nx_;
   int ny = sim.ny_;
   int nz = sim.nz_;
   
   vector<Long64> myCells(anatomy.size());
   for (unsigned ii=0; ii<anatomy.size(); ++ii)
      myCells[ii] = anatomy.gid(ii);
   
   GridRouter router(myCells, nx, ny, nz, MPI_COMM_WORLD);
   sim.sendMap_ = router.sendMap();
   sim.commTable_ = new CommTable(router.commTable());

   HaloExchange<AnatomyCell> cellExchange(sim.sendMap_, (sim.commTable_));

   int recvBuffSize=cellExchange.recvSize(); 
   int nRemote = recvBuffSize;
   int nLocal = anatomy.nLocal();
   {
      lazy_array<AnatomyCell> lzCells;
      lzCells.resize(nLocal+nRemote);
      wo_array_ptr<AnatomyCell> wo_ptr = lzCells.writeonly(CPU);
      copy(anatomy.cellArray().begin(), anatomy.cellArray().end(), wo_ptr.begin());
      cellExchange.execute(lzCells, nLocal);

      anatomy.cellArray().resize(lzCells.size());
      ro_array_ptr<AnatomyCell> ro_ptr = lzCells.readonly(CPU);
      copy(ro_ptr.begin(), ro_ptr.end(), anatomy.cellArray().begin());
   }
   anatomy.nRemote() = nRemote;

   OBJECT* obj = object_find(name.c_str(), "DECOMPOSITION");
   bool printStats;
   objectGet(obj, "printStats", printStats, "0");
   if (printStats)
   {
      int myRank; 
      MPI_Comm_rank(MPI_COMM_WORLD, & myRank);
      PFILE* file = Popen("nNbrTasks", "w", comm);
      Pprintf(file, "%d : %d %d %d %d\n", 
      myRank,anatomy.nGlobal(),anatomy.nLocal(),anatomy.nRemote(),sim.commTable_->nNbrs());
      Pclose(file);

      file = Popen("msgSizes", "w", comm);
      vector<int> mSize = sim.commTable_->msgSize();
      for (unsigned ii=0; ii<mSize.size(); ++ii) Pprintf(file, "%d\n", mSize[ii]);
      Pclose(file);
   }
   

}

