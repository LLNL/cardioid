#include "AnatomyReader.hh"

#include <cassert>
#include <set>
#include <vector>
#include "pio.h"
#include "object_cc.hh"
#include "Anatomy.hh"
#include "BucketOfBits.hh"
#include "readPioFile.hh"

using std::string;
using std::set;
using std::vector;

namespace
{
   void readBucket(Anatomy& anatomy, const set<int>& typeSet,
                   BucketOfBits* bucketP);
}

BucketOfBits* readAnatomy(const string& filename, MPI_Comm comm,
                          Anatomy& anatomy, const set<int>& typeSet)
{
   PFILE* file = Popen(filename.c_str(), "r", comm);

   OBJECT* hObj = file->headerObject;

   int nx, ny, nz;
   objectGet(hObj, "nx", nx, "0");
   objectGet(hObj, "ny", ny, "0");
   objectGet(hObj, "nz", nz, "0");
//   objectGet(hObj, "cellTypes", cellType_s);

   Long64 nGlobal = nx*ny*nz;
   assert(nGlobal > 0);
//   assert(cellType_s.size() > 0);

   anatomy.setGridSize(nx, ny,nz);

   BucketOfBits* bucketP = readPioFile(file);
   Pclose(file);

   readBucket(anatomy, typeSet, bucketP);
   return bucketP;
}

namespace
{
   void readBucket(Anatomy& anatomy, const set<int>& typeSet,
                   BucketOfBits* bucketP)
   {
      vector<AnatomyCell>& cells = anatomy.cellArray();
      
      unsigned nRecords =  bucketP->nRecords();
      unsigned nFields =   bucketP->nFields();
      unsigned gidIndex =  bucketP->getIndex("gid");
      unsigned typeIndex = bucketP->getIndex("cellType");
      assert( gidIndex != nFields &&
              typeIndex != nFields);
      for (unsigned ii=0; ii<nRecords; ++ii)
      {
         AnatomyCell tmp;
         BucketOfBits::Record rr = bucketP->getRecord(ii);
         rr.getValue(gidIndex, tmp.gid_);
         rr.getValue(typeIndex, tmp.cellType_);
         if (typeSet.count(tmp.cellType_) != 0)
            cells.push_back(tmp);
      }
   }
}
