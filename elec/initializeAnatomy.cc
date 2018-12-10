#include "initializeAnatomy.hh"

#include <iostream>
#include <cassert>

#include "AnatomyReader.hh"
#include "object_cc.hh"
#include "Anatomy.hh"
#include "TupleToIndex.hh"
#include "setConductivity.hh"
#include "BucketOfBits.hh"
#include "Prand48Object.hh"

using namespace std;

namespace
{
   // Caller is resposible to delete the returned pointer
   BucketOfBits* readUsingPio(Anatomy& anatomy,
                              OBJECT* obj, MPI_Comm comm);
   BucketOfBits* generateTissueBrick(Anatomy& anatomy, OBJECT* obj, MPI_Comm comm);
}

/*!
  @page obj_ANATOMY ANATOMY object
  Defines the anatomical structure to simulate.

  @beginkeywords
    @kw{conductivity, The name of the CONDUCTIVITY object which specifies
      the method and parameters to determine cell conductivity., No Default}
    @kw{dx, Cell size in the x-direction., 0.2 mm}
    @kw{dy, Cell size in the y-direction., 0.2 mm}
    @kw{dz, Cell size in the z-direction., 0.2 mm}
    @kw{method, Choose from "brick" or "pio", pio}
    @endkeywords

  @subpage ANATOMY_brick

  @subpage ANATOMY_pio

*/
void initializeAnatomy(Anatomy& anatomy, const string& name, MPI_Comm comm)
{
   OBJECT* obj = object_find(name.c_str(), "ANATOMY");

   string conductivityName;
   objectGet(obj, "conductivity", conductivityName, "");
   double dx, dy, dz;
   objectGet(obj, "dx", dx, "0.2", "l");
   objectGet(obj, "dy", dy, "0.2", "l");
   objectGet(obj, "dz", dz, "0.2", "l");
   anatomy.dx() = dx;
   anatomy.dy() = dy;
   anatomy.dz() = dz;
   double offset_x, offset_y, offset_z;
   objectGet(obj, "offset_x", offset_x, "0.0", "l");
   objectGet(obj, "offset_y", offset_y, "0.0", "l");
   objectGet(obj, "offset_z", offset_z, "0.0", "l");
   anatomy.offset_x() = offset_x;
   anatomy.offset_y() = offset_y;
   anatomy.offset_z() = offset_z;


   BucketOfBits* data = 0;

   string method;
   objectGet(obj, "method", method, "pio");
   if (method == "pio")
      data = readUsingPio(anatomy, obj, comm);
   else if (method == "brick")
      data = generateTissueBrick(anatomy, obj, comm);
   else if (method == "simple")
      // We can wire in the simple load code that Erik originally wrote
      // here if we still need it.
      assert(1==0);
   else
      assert(1==0);

   Long64 nGlobal;
   Long64 nLocal = anatomy.size();
   MPI_Allreduce(&nLocal, &nGlobal, 1, MPI_LONG_LONG, MPI_SUM, comm);
   anatomy.nGlobal() = nGlobal;

   Tuple globalGridSize(anatomy.nx(), anatomy.ny(), anatomy.nz());

   //ewd: check that there are no gids greater than or equal to nx*ny*nz, which
   //ewd: will cause the load balancers to hang or crash.
   Long64 nXYZ = anatomy.nx()*anatomy.ny()*anatomy.nz();
   for (unsigned ii=0; ii<nLocal; ++ii)
   {
      Long64 gid = anatomy.gid(ii);
      if (gid >= nXYZ || gid < 0)
      {
         cout << "Invalid gid " << gid << " for grid " << anatomy.nx() << " x " << anatomy.nx() << " x " << anatomy.nx() << " (" << nXYZ << " points total)" << endl;
         assert(gid >= 0);
         assert(gid < nXYZ);
      }
   }

   setConductivity(conductivityName, *data, globalGridSize, anatomy.cellArray());
   delete data;
}






namespace
{
   /*!
     @page ANATOMY_pio ANATOMY pio method

     Specifies that the anatomy should be read from a file.  The file
     format for the anatomy file is described in
     section @ref sec_anatomy_format.

     @issue{Why don't we store the cell size (or at least the default cell
     size) in the anatomy file?}

     @beginkeywords
     @kw{filename, Path to pio file to read.  Includes the #-sign\, but not
        a six digit sequence number., snapshot.initial/anatomy#}
     @endkeywords
   */
   BucketOfBits* readUsingPio(Anatomy& anatomy,
                              OBJECT* obj, MPI_Comm comm)
   {
      int myRank;
      MPI_Comm_rank(comm, &myRank);

      string fileName;
      objectGet(obj, "fileName", fileName, "snapshot.initial/anatomy#");

      if (myRank==0) cout << "Starting read" <<endl;

      BucketOfBits* bucketP = readAnatomy(fileName, comm, anatomy);
      if (myRank==0) cout << "Finished read" <<endl;
      return bucketP;
   }
}


namespace
{
   /*!
     @page ANATOMY_brick ANATOMY brick method

     Generates an orthorhombic "brick" of cells of a specified size.
     All cells are of the same type.

     @issue{Something should be said about the layout of the
     discretization with respect to the specfied box.}

     @beginkeywords
       @kw{cellType, The cell type in the brick, 102}
       @kw{xSize, Size of the simulation in the x-direction in mm, 3 mm}
       @kw{ySize, Size of the simulation in the y-direction in mm, 7 mm}
       @kw{zSize, Size of the simulation in the z-direction in mm, 20 mm}
     @endkeywords

     ~~~~
     niederer ANATOMY
     {
       method = brick;
       cellType = 102;
       dx = 0.1;
       dy = 0.1;
       dz = 0.1;
       xSize = 3;
       ySize = 7;
       zSize = 20;
     }
     ~~~~

   */
   BucketOfBits*  generateTissueBrick(Anatomy& anatomy, OBJECT* obj, MPI_Comm comm)
   {
      int myRank;
      int nTasks;
      MPI_Comm_size(comm, &nTasks);
      MPI_Comm_rank(comm, &myRank);

      bool sizeSpecified = (object_testforkeyword(obj, "xSize") ||
                            object_testforkeyword(obj, "ySize") ||
                            object_testforkeyword(obj, "zSize") );
      bool numberSpecified = (object_testforkeyword(obj, "nx") ||
                              object_testforkeyword(obj, "ny") ||
                              object_testforkeyword(obj, "nz") );
      double xSize, ySize, zSize;
      objectGet(obj, "xSize", xSize, "3", "l");
      objectGet(obj, "ySize", ySize, "7", "l");
      objectGet(obj, "zSize", zSize, "20","l");

      int cellType;
      string cellString;
      uint64_t seed;
      objectGet(obj, "seed",     seed,       "195782659275");
      objectGet(obj, "cellType", cellString, "102");
      if (cellString != "random")
         cellType = atoi(cellString.c_str());

      int nx = int(xSize/anatomy.dx());
      int ny = int(ySize/anatomy.dy());
      int nz = int(zSize/anatomy.dz());

      if (numberSpecified)
      {
         assert(!sizeSpecified);
         objectGet(obj, "nx", nx, "16");
         objectGet(obj, "ny", ny, "16");
         objectGet(obj, "nz", nz, "14");
      }

      anatomy.setGridSize(nx, ny, nz);

      Long64 maxGid = Long64(nx)*Long64(ny)*Long64(nz);
      unsigned cellsPerTask = maxGid/nTasks;
      if (maxGid%nTasks != 0)
         ++cellsPerTask;

      Long64 gidBegin = min(maxGid, Long64(myRank)*Long64(cellsPerTask));
      Long64 gidEnd   = min(maxGid, gidBegin+Long64(cellsPerTask));
      unsigned nLocal = gidEnd - gidBegin;

      vector<AnatomyCell>& cells = anatomy.cellArray();
      cells.reserve(nLocal);
      for (Long64 ii=gidBegin; ii<gidEnd; ++ii)
      {
         AnatomyCell tmp;
         tmp.gid_ = ii;
         if (cellString == "random")
         {
            Prand48Object rand(tmp.gid_, seed, 0xace2468bdf1357llu);
            tmp.cellType_ = 100 + rand(3);
         }
         else
            tmp.cellType_ = cellType;
         cells.push_back(tmp);
      }
      return new BucketOfBits(vector<string>(), vector<string>(), vector<string>());
   }
}
