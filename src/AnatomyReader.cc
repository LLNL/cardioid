#include "AnatomyReader.hh"

#include <cassert>
#include "pio.h"
#include "pioFixedRecordHelper.h"
#include "object_cc.hh"

using std::string;


AnatomyReader::AnatomyReader(const string& filename, MPI_Comm comm)
{
   PFILE* file = Popen(filename.c_str(), "r", comm);

   OBJECT* hObj = file->headerObject;

   objectGet(hObj, "nx", _nx, "0");
   objectGet(hObj, "ny", _ny, "0");
   objectGet(hObj, "nz", _nz, "0");
   objectGet(hObj, "cellTypes", _cellTypes);

   assert(_nx*_ny*_nz > 0);
//   assert(_cellTypes.size() > 0);
   
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
   
   Pclose(file);
}

void AnatomyReader::asciiReader(PFILE* file)
{
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
      _anatomy.push_back(tmp);
   }
}

void AnatomyReader::binaryReader(PFILE* file)
{
   assert(1==0);
}
