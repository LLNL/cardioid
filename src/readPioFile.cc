#include "readPioFile.hh"

#include <cassert>
#include "pio.h"
#include "pioFixedRecordHelper.h"
#include "pioVariableRecordHelper.h"
#include "object_cc.hh"
#include "BucketOfBits.hh"
#include "ioUtils.h"
#include <cstring>

using namespace std;


namespace
{
   void fixRecordAscii (PFILE* file, BucketOfBits* bucketP);
   void fixRecordBinary(PFILE* file, BucketOfBits* bucketP);
   void varRecordAscii (PFILE* file, BucketOfBits* bucketP);
   void readAscii (PFILE* file, unsigned nRecords, BucketOfBits* bucketP);
   void readBinary(PFILE* file, unsigned lrec, unsigned nRecords, BucketOfBits* bucketP);
}



/** We can always count on a pio file to have field_names and
 *  field_types, but field_units are optional.  If there are no units
 *  then we just claim that all fields are dimensionless.
 *
 *  If we read an ascii file, the stl string that represents the record
 *  in the BucketOfBits will include the \n that ends the ascii record.
 */
BucketOfBits* readPioFile(PFILE* file)
{
   OBJECT* hObj = file->headerObject;
   vector<string> fieldNames;
   vector<string> fieldTypes;
   vector<string> fieldUnits;
   objectGet(hObj, "field_names", fieldNames);
   objectGet(hObj, "field_types", fieldTypes);
   if (object_testforkeyword(hObj, "field_units"))
      objectGet(hObj, "field_units", fieldUnits);
   else
      fieldUnits.assign(fieldNames.size(), "1");
   
   BucketOfBits* bucketP = new BucketOfBits(fieldNames, fieldTypes, fieldUnits);

   switch (file->datatype)
   {
     case FIXRECORDASCII:
      fixRecordAscii(file, bucketP);
      break;
     case FIXRECORDBINARY:
      fixRecordBinary(file, bucketP);
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
   /** Handling of the endian key is not completely satisfactory as it
    * is not thread safe.  If we ever happen to have two different
    * threads reading two different binary files with different
    * endian-ness we will have a problem.  */
   void fixRecordBinary(PFILE* file, BucketOfBits* bucketP)
   {
      PIO_FIXED_RECORD_HELPER* helper = (PIO_FIXED_RECORD_HELPER*) file->helper;
      unsigned lrec = helper->lrec;
      unsigned nRecords = file->bufsize/lrec;
      assert(file->bufsize%lrec == 0);
      OBJECT* hObj = file->headerObject;
      unsigned key;
      objectGet(hObj, "endian_key", key, "0");
      assert(key != 0); // This can only fail if key isn't set in file.
      ioUtils_setSwap(key);
      readBinary(file, lrec, nRecords, bucketP);
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

namespace
{
   void readBinary(PFILE* file, unsigned lrec, unsigned nRecords, BucketOfBits* bucketP)
   {
      char buf[lrec];
      for (unsigned ii=0; ii<nRecords; ++ii)
      {
         Pread(buf, lrec, 1, file);
         bucketP->addRecord(string(buf, lrec));
      }
   }
}
