// $Id$

#ifndef PIO_VARIABLE_RECORD_HELPER_H
#define PIO_VARIABLE_RECORD_HELPER_H

#include "pioHelper.h"

#ifdef __cplusplus
extern "C" {
#endif

/** ASCII Version */
typedef struct PioVariableRecordAsciiHelper_st
{
   phb_destroy destroy;
   phb_endOfRecords endOfRecords;
   // Items for "inherited" class
   char delimiter;
} PIO_VARIABLE_RECORD_ASCII_HELPER;

PIO_HELPER* pvrah_create(OBJECT* header);
unsigned pvrah_nRecords(const char* buf, size_t bufsize, const char delim);


/** BINARY version */
typedef struct PioVariableRecordBinaryHelper_st
{
   phb_destroy destroy;
   phb_endOfRecords endOfRecords;
   // Items for "inherited" class
   char* lrecType; // record type of lrec field (i.e., b1)
   unsigned lrecSize; // number of bytes needed to specify record length
   unsigned offset; // offset from start of record to lrec in bytes.
} PIO_VARIABLE_RECORD_BINARY_HELPER;

PIO_HELPER* pvrbh_create(OBJECT* header);

#ifdef __cplusplus
}
#endif
#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
