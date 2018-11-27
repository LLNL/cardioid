// $Id$

#ifndef PIO_HELPER_H
#define PIO_HELPER_H

#include "object.h"

#ifdef __cplusplus
extern "C" {
#endif

struct PioHelper_st;


typedef void    (*phb_destroy)      (struct PioHelper_st* pioHelper);
typedef size_t  (*phb_endOfRecords) (const char* buf,
				     size_t bufSize,
				     struct PioHelper_st* pioHelper);

/** If this was C++ this would be an abstract base class.  The purpose
 *  of the PIO_HELPER struct is to define an interface that a pfile of a
 *  particular datatype must support to be readable by pio.  Note that
 *  datatype is required to be defined in the pfile header.  For an
 *  example of a concrete class see pioFixedRecordHelper.c.  Each
 *  concrete class must provide an implementation of the following
 *  functions:
 *
 *  destroy:  This is the virtual destructor.  Anyone who creates
 *    a PIO_HELPER is expected to call destroy when the helper is
 *    no longer needed.  This function should free any internal
 *    storage needed by the concrete helper class.
 *
 *  endOfRecords: When pio reads a data from disk it is required to
 *    send an integer number of records to each tasks in the I/O
 *    group.  However, for an arbitrary file type pio has no way of
 *    knowing where a record end is.  To resolve this issue, each
 *    time a chunk of data is read from the disk
 *    the reader task pass the buffer to endOfRecords.  This
 *    function will analyze the buffer and return the location of
 *    the end of a complete record in the buffer.  Pio will send
 *    the complete records to a task in the I/O group and prepend
 *    any leftover bytes to the start of the next chunk read from
 *    disk.  This way each task is guaranteed to receive an integer
 *    number of records.
 *    
 *  create: This function isn't in in the abstract base class,
 *    but would be if this were C++ and we had constructors.
 *    The pioHelperFactory will call the create function to
 *    initialize the internals of any given concrete class.
 *    Since the only information passed to the factory is the
 *    file header (as an OBJECT*) the only practical signature
 *    for the create function is create(OBJECT* header);
 *
 *  Each concrete class must also provide appropriate implementation in
 *  the pioHelperFactory function to create the concrete class when the
 *  corresponding value of datatype is found in the pfile header.
*/
typedef struct PioHelper_st
{
   phb_destroy destroy;
   phb_endOfRecords endOfRecords;
} PIO_HELPER;

/** This factory function creates a concrete helper class that
 *  corresponds to the datatype specified in the header. */
PIO_HELPER* pioHelperFactory(OBJECT* header);


#ifdef __cplusplus
}
#endif

#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
