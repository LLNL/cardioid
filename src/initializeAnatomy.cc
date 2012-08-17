#include "initializeAnatomy.hh"

#include <iostream>
#include <cassert>

#include "AnatomyReader.hh"
#include "object_cc.hh"
#include "Anatomy.hh"
#include "TupleToIndex.hh"
#include "setConductivity.hh"
#include "BucketOfBits.hh"
#include "Drand48Object.hh"

using namespace std;

namespace
{
   // Caller is resposible to delete the returned pointer
   BucketOfBits* readUsingPio(Anatomy& anatomy, const set<int>& typeSet,
                              OBJECT* obj, MPI_Comm comm);
   BucketOfBits* generateTissueBrick(Anatomy& anatomy, OBJECT* obj, MPI_Comm comm);
}

   

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

   set<int> typeSet;
   typeSet.insert(100);
   typeSet.insert(101);
   typeSet.insert(102);
   typeSet.insert(30);
   typeSet.insert(31);
   typeSet.insert(32);
   typeSet.insert(33);
   typeSet.insert(34);
   typeSet.insert(35);
   typeSet.insert(75);
   typeSet.insert(76);
   typeSet.insert(77);
   
   BucketOfBits* data = 0;
   
   string method;
   objectGet(obj, "method", method, "pio");
   if (method == "pio")
      data = readUsingPio(anatomy, typeSet, obj, comm);
   else if (method == "brick")
      data = generateTissueBrick(anatomy, obj, comm);
   else if (method == "simple")
      // We can wire in the simple load code that Erik originally wrote
      // here if we still need it.
      assert(1==0);
   else
      assert(1==0);

   int nGlobal;
   int nLocal = anatomy.size();
   MPI_Allreduce(&nLocal, &nGlobal, 1, MPI_INT, MPI_SUM, comm);
   anatomy.nGlobal() = nGlobal;
   
   Tuple globalGridSize(anatomy.nx(), anatomy.ny(), anatomy.nz());
   
   setConductivity(conductivityName, *data, globalGridSize, typeSet, anatomy.cellArray());
   delete data;
}


namespace
{
   BucketOfBits* readUsingPio(Anatomy& anatomy, const set<int>& typeSet,
                              OBJECT* obj, MPI_Comm comm)
   {
      int myRank;
      MPI_Comm_rank(comm, &myRank);
      
      string fileName;
      objectGet(obj, "fileName", fileName, "snapshot.initial/anatomy#");
      
      if (myRank==0) cout << "Starting read" <<endl;

      BucketOfBits* bucketP = readAnatomy(fileName, comm, anatomy, typeSet);
      if (myRank==0) cout << "Finished read" <<endl;
      return bucketP;
   }
}

namespace
{
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
      objectGet(obj, "cellType", cellString, "102");
      if (cellString != "random")
         cellType = atoi(cellString.c_str());
      Drand48Object rand(myRank);  
      
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
            tmp.cellType_ = 100 + rand(3);
         else
            tmp.cellType_ = cellType;
         cells.push_back(tmp);
      }
      return new BucketOfBits(vector<string>(), vector<string>(), vector<string>());
   }
}

