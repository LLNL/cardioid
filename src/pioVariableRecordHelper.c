#include "pioVariableRecordHelper.h"

#include <assert.h>
#include <string.h>
#include "ddcMalloc.h"
#include "ioUtils.h"

// ASCII
static void    pvrah_destroy(PIO_HELPER* this);
static size_t  pvrah_endOfRecords(const char* buf,
				  size_t nBuf,
				  PIO_HELPER* thisBase);

// BINARY
static void    pvrbh_destroy(PIO_HELPER* this);
static size_t  pvrbh_endOfRecords(const char* buf,
				  size_t nBuf,
				  PIO_HELPER* thisBase);

// ASCII
PIO_HELPER* pvrah_create(OBJECT* header)
{
   PIO_VARIABLE_RECORD_ASCII_HELPER* helper =
      ddcMalloc(sizeof(PIO_VARIABLE_RECORD_ASCII_HELPER));

   helper->destroy = pvrah_destroy;
   helper->endOfRecords = pvrah_endOfRecords;
   helper->delimiter = '\n';

   return (PIO_HELPER*) helper;
}

void pvrah_destroy(PIO_HELPER* this)
{
   return;
}

size_t pvrah_endOfRecords(const char* buf, size_t nBuf, PIO_HELPER* thisBase)
{
   PIO_VARIABLE_RECORD_ASCII_HELPER* this = (PIO_VARIABLE_RECORD_ASCII_HELPER*) thisBase;
   
   size_t nBytes=nBuf;
   while (nBytes > 0)
   {
      if (buf[nBytes-1] == this->delimiter) break;
      --nBytes;
   }
   return nBytes;
}

unsigned pvrah_nRecords(const char* buf, size_t bufsize, const char delim)
{
   unsigned nRecords = 0;
   for (size_t ii=0; ii<bufsize; ++ii) if (buf[ii] == delim) ++nRecords;

   return nRecords;
}


// BINARY

/** To support varible length binary records, there must be a field in
 *  the record that gives the record length.  It does not need to be the
 *  first field (for example, the checksum may come first), but it must
 *  be named lrec, and the record length should include the entire
 *  record, including any fields that come before the lrec field.
 *
 *  This function identifies which field contains the lrec data, the
 *  type of that field, and the offset in bytes from the start of the
 *  record.  This is all the information that is needed by the
 *  endOfRecords function to scan through data to find the last valid
 *  record.
 */
PIO_HELPER* pvrbh_create(OBJECT* header)
{
   PIO_VARIABLE_RECORD_BINARY_HELPER* helper =
      ddcMalloc(sizeof(PIO_VARIABLE_RECORD_BINARY_HELPER));

   helper->destroy = pvrbh_destroy;
   helper->endOfRecords = pvrbh_endOfRecords;

   unsigned endianKey;
   object_get(header, "endian_key", &endianKey, INT, 1, "0");
   ioUtils_setSwap(endianKey);
   
   char** fieldNames = NULL;
   char** fieldTypes = NULL;
   unsigned nFields =
   object_getv(header, "field_names", (void*)&fieldNames, STRING,ABORT_IF_NOT_FOUND);
   object_getv(header, "field_types", (void*)&fieldTypes, STRING,ABORT_IF_NOT_FOUND);

   int lrecField = -1;
   for (unsigned ii=0; ii<nFields; ++ii)
   if (strcmp(fieldNames[ii], "lrec") == 0)
   {
      lrecField = ii;
      break;
   }

   assert(lrecField >= 0);
   
   helper->lrecSize = atoi(fieldTypes[lrecField]+1);
   helper->lrecType = strdup(fieldTypes[lrecField]);
   helper->offset = 0;
   for (int ii=0; ii<lrecField; ++ii)
      helper->offset += atoi(fieldTypes[ii]+1);

   for (unsigned ii=0; ii<nFields; ++ii)
   {
      ddcFree(fieldTypes[ii]);
      ddcFree(fieldNames[ii]);
   }
   ddcFree(fieldNames);
   ddcFree(fieldTypes);
   
   return (PIO_HELPER*) helper;
}

void pvrbh_destroy(PIO_HELPER* thisBase)
{
   PIO_VARIABLE_RECORD_BINARY_HELPER* this =
      (PIO_VARIABLE_RECORD_BINARY_HELPER*) thisBase;
   
   ddcFree(this->lrecType);
   return;
}

/** This routine assumes that the first byte of buf is the start of a
 *  record.  It finds the length each record and jumps ahead by that
 *  amount until only a partial record remains.  */
size_t pvrbh_endOfRecords(const char* buf,
			 size_t nBuf,
			 PIO_HELPER* thisBase)
{
   PIO_VARIABLE_RECORD_BINARY_HELPER* this =
      (PIO_VARIABLE_RECORD_BINARY_HELPER*) thisBase;
   
   size_t eor = 0;

   // The conditional for this while loop catches the case where there
   // the record length field isn't present within the remaining bytes.
   while (eor + this->offset + this->lrecSize < nBuf)
   {
      size_t lrec = mkInt((const unsigned char*) buf+eor+this->offset, this->lrecType);
      if (eor + lrec > nBuf)
	 break;
      eor += lrec;
   }

   return eor;
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
