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
   void generateTissueBrick(Simulate& sim, OBJECT* obj);
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
      generateTissueBrick(sim, obj);
   else if (method == "simple")
      // We can wire in the simple load code that Erik originally wrote
      // here if we still need it.
      assert(1==0);
   else
      assert(1==0);

   sim.anatomy_.setGridSize(sim.nx_, sim.ny_, sim.nz_);
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
   void generateTissueBrick(Simulate& sim, OBJECT* obj)
   {
      double xSize, ySize, zSize;
      int cellType;
      objectGet(obj, "xSize", xSize, "3");
      objectGet(obj, "ySize", ySize, "7");
      objectGet(obj, "zSize", zSize, "20");
      objectGet(obj, "cellType", cellType, "102");

      Anatomy& anatomy = sim.anatomy_;
      
      int nx = xSize/anatomy.dx();
      int ny = ySize/anatomy.dy();
      int nz = zSize/anatomy.dz();
      anatomy.setGridSize(nx, ny, nz);

      unsigned nCells = nx*ny*nz;
      vector<AnatomyCell> cells = sim.anatomy_.cellArray();
      cells.reserve(nCells);
      TupleToIndex tupleToIndex(nx, ny, nz);
      
      for (unsigned iz=0; iz<nz; ++iz)
         for (unsigned iy=0; iy<ny; ++iy)
            for (unsigned ix=0; ix<nx; ++ix)
            {
               AnatomyCell tmp;
               tmp.gid_ = tupleToIndex(ix, iy, iz);
               tmp.theta_= 0.0;
               tmp.phi_ = 0.0;
               tmp.cellType_ = cellType;
               cells.push_back(tmp);
            }
   }
}
