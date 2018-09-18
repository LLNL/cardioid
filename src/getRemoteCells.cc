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
   anatomy.nRemote() = recvBuffSize;
   anatomy.cellArray().resize(anatomy.cellArray().size()+anatomy.nRemote());
   //HACK!  We have a cpu-only data structure, but that should be fine....
   //Make the data fit the interface.  
   rw_larray_ptr<AnatomyCell> CPUonlyCellArray(NULL, &(anatomy.cellArray()[0]), 0, anatomy.cellArray().size());
   cellExchange.execute(CPUonlyCellArray, anatomy.nLocal());
   //size_t oldSize = anatomy.size();
   //cellExchange.execute(anatomy.cellArray(), anatomy.nLocal());
   //anatomy.nRemote() = anatomy.size() - oldSize;
//   dumpCells(anatomy);

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

