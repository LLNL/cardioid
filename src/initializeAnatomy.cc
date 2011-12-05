#include "initializeAnatomy.hh"

#include <iostream>
#include <cassert>

#include "Simulate.hh"
#include "AnatomyReader.hh"
#include "object_cc.hh"
#include "Anatomy.hh"
#include "TupleToIndex.hh"
using namespace std;

namespace
{
   void readUsingPio(Simulate& sim, OBJECT* obj, MPI_Comm comm);
   void generateTissueBrick(Simulate& sim, OBJECT* obj, MPI_Comm comm);
}

   

void initializeAnatomy(Simulate& sim, const string& name, MPI_Comm comm)
{
   OBJECT* obj = object_find(name.c_str(), "ANATOMY");

   double dx, dy, dz;
   objectGet(obj, "dx", dx, "0.2");
   objectGet(obj, "dy", dy, "0.2");
   objectGet(obj, "dz", dz, "0.2");
   sim.anatomy_.dx() = dx;
   sim.anatomy_.dy() = dy;
   sim.anatomy_.dz() = dz;

   string method;
   objectGet(obj, "method", method, "pio");

   if (method == "pio")
      readUsingPio(sim, obj, comm);
   else if (method == "brick")
      generateTissueBrick(sim, obj, comm);
   else if (method == "simple")
      // We can wire in the simple load code that Erik originally wrote
      // here if we still need it.
      assert(1==0);
   else
      assert(1==0);
}


namespace
{
   void readUsingPio(Simulate& sim, OBJECT* obj, MPI_Comm comm)
   {
      int myRank;
      MPI_Comm_rank(comm, &myRank);
      
      string fileName;
      objectGet(obj, "fileName", fileName, "snapshot.initial/atatomy#");
      
      if (myRank==0) cout << "Starting read" <<endl;

      AnatomyReader reader(fileName, comm, sim);
      if (myRank==0) cout << "Finished read" <<endl;
   }
}

namespace
{
   void generateTissueBrick(Simulate& sim, OBJECT* obj, MPI_Comm comm)
   {
      int myRank;
      int nTasks;
      MPI_Comm_size(comm, &nTasks);
      MPI_Comm_rank(comm, &myRank);

      double xSize, ySize, zSize;
      int cellType;
      objectGet(obj, "xSize", xSize, "3");
      objectGet(obj, "ySize", ySize, "7");
      objectGet(obj, "zSize", zSize, "20");
      objectGet(obj, "cellType", cellType, "102");

      Anatomy& anatomy = sim.anatomy_;
      
      int nx = int(xSize/anatomy.dx());
      int ny = int(ySize/anatomy.dy());
      int nz = int(zSize/anatomy.dz());
      anatomy.setGridSize(nx, ny, nz);

      Long64 maxGid = Long64(nx)*Long64(ny)*Long64(nz);
      unsigned cellsPerTask = maxGid/nTasks;
      if (maxGid%nTasks != 0)
         ++cellsPerTask;

      Long64 gidBegin = min(maxGid, Long64(myRank)*Long64(cellsPerTask));
      Long64 gidEnd   = min(maxGid, gidBegin+Long64(cellsPerTask));
      unsigned nLocal = gidEnd - gidBegin;

      vector<AnatomyCell>& cells = sim.anatomy_.cellArray();
      cells.reserve(nLocal);
      for (Long64 ii=gidBegin; ii<gidEnd; ++ii)
      {
         AnatomyCell tmp;
         tmp.gid_ = ii;
         tmp.theta_= 0;
         tmp.phi_ = 0;
         tmp.cellType_ = cellType;
         cells.push_back(tmp);
      }
   }
}
