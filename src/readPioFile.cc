#include "readPioFile.hh"

#include <cassert>
#include "pio.h"
#include "pioFixedRecordHelper.h"
#include "pioVariableRecordHelper.h"
#include "object_cc.hh"
#include "BucketOfBits.hh"

using namespace std;


namespace
{
   void fixRecordAscii(PFILE* file, BucketOfBits* bucketP);
   void varRecordAscii(PFILE* file, BucketOfBits* bucketP);
   void readAscii(PFILE* file, unsigned nRecords, BucketOfBits* bucketP);
}



BucketOfBits* readPioFile(PFILE* file)
{
   OBJECT* hObj = file->headerObject;
   vector<string> fieldNames;
   vector<string> fieldTypes;
   objectGet(hObj, "field_names", fieldNames);
   objectGet(hObj, "field_types", fieldTypes);

   BucketOfBits* bucketP = new BucketOfBits(fieldNames, fieldTypes);

   switch (file->datatype)
   {
     case FIXRECORDASCII:
      fixRecordAscii(file, bucketP);
      break;
     case VARRECORDASCII:
      varRecordAscii(file, bucketP);
      break;
     default:
      assert(false);
   }

   return bucketP;
}

namespace
{
   void fixRecordAscii(PFILE* file, BucketOfBits* bucketP)
   {
      PIO_FIXED_RECORD_HELPER* helper = (PIO_FIXED_RECORD_HELPER*) file->helper;
      unsigned lrec = helper->lrec;
      unsigned nRecords = file->bufsize/lrec;
      assert(file->bufsize%lrec == 0);
      readAscii(file, nRecords, bucketP);
   }
}

namespace
{
   void varRecordAscii(PFILE* file, BucketOfBits* bucketP)
   {
      PIO_VARIABLE_RECORD_ASCII_HELPER* helper =
         (PIO_VARIABLE_RECORD_ASCII_HELPER*) file->helper;
      unsigned nRecords = pvrah_nRecords(file->buf, file->bufsize, helper->delimiter);
      readAscii(file, nRecords, bucketP);
   }
}

namespace
{
   void readAscii(PFILE* file, unsigned nRecords, BucketOfBits* bucketP)
   {
      unsigned maxRec = 2048;
      char buf[maxRec];
      for (unsigned ii=0; ii<nRecords; ++ii)
      {
         Pfgets(buf, maxRec, file);
         assert(strlen(buf) < maxRec);
         bucketP->addRecord(buf);
      }
   }
}
