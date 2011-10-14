#include "AnatomyReader.hh"

#include <cassert>
#include <set>
#include "pio.h"
#include "pioFixedRecordHelper.h"
#include "object_cc.hh"
#include "Simulate.hh"

using std::string;
using std::set;


AnatomyReader::AnatomyReader(const string& filename, MPI_Comm comm,
			     Simulate& sim)
: _anatomy(sim.cell_)
{
   int myRank;
   MPI_Comm_rank(comm, &myRank);

   PFILE* file = Popen(filename.c_str(), "r", comm);

   OBJECT* hObj = file->headerObject;

   int nx, ny, nz;
   objectGet(hObj, "nx", nx, "0");
   objectGet(hObj, "ny", ny, "0");
   objectGet(hObj, "nz", nz, "0");
//   objectGet(hObj, "cellTypes", _cellTypes);

   assert(nx*ny*nz > 0);
//   assert(_cellTypes.size() > 0);

   sim.nx_ = nx;
   sim.ny_ = ny;
   sim.nz_ = nz;
   
   
   switch (file->datatype)
   {
     case FIXRECORDASCII:
      asciiReader(file);
      break;
     case FIXRECORDBINARY:
      binaryReader(file);
      break;
     default:
      assert(1==0);
   }

   for (unsigned ii=0; ii<_anatomy.size(); ++ii)
     _anatomy[ii]._dest = myRank;
   
   Pclose(file);
}

void AnatomyReader::asciiReader(PFILE* file)
{
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
   

   PIO_FIXED_RECORD_HELPER* helper = (PIO_FIXED_RECORD_HELPER*) file->helper;
   unsigned lrec = helper->lrec;
   unsigned nRecords = file->bufsize/lrec;
   assert(file->bufsize%lrec == 0);
   
   for (unsigned ii=0; ii<nRecords; ++ii)
   {
      char buf[file->recordLength];
      AnatomyCell tmp;
      buf[0] = '\0';
      Pfgets(buf, file->recordLength, file);
      sscanf(buf, "%llu %d %d %d", &(tmp._gid), &(tmp._cellType), &(tmp._theta), &(tmp._phi));
      if (typeSet.find(tmp._cellType) != typeSet.end())
	 _anatomy.push_back(tmp);
   }
}

void AnatomyReader::binaryReader(PFILE* file)
{
   assert(1==0);
}
