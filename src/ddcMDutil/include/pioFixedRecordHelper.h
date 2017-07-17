// $Id$

#ifndef PIO_FIXED_RECORD_HELPER_H
#define PIO_FIXED_RECORD_HELPER_H

#include "pioHelper.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct PioFixedRecordHelper_st
{
   phb_destroy destroy;
   phb_endOfRecords endOfRecords;
   // Items for "inherited" class
   unsigned lrec;
} PIO_FIXED_RECORD_HELPER;

PIO_HELPER* pfrh_create(OBJECT* header);

#ifdef __cplusplus
}
#endif
#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
