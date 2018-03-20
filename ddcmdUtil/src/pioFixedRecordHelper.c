#include "pioFixedRecordHelper.h"

#include <assert.h>
#include "ddcMalloc.h"

static void    pfrh_destroy(PIO_HELPER* this);
static size_t  pfrh_endOfRecords(const char* buf,
				 size_t nBuf,
				 PIO_HELPER* this);


PIO_HELPER* pfrh_create(OBJECT* header)
{
   PIO_FIXED_RECORD_HELPER* helper =
      ddcMalloc(sizeof(PIO_FIXED_RECORD_HELPER));

   object_get(header, "lrec", &helper->lrec, INT, 1, "0");
   assert(helper->lrec > 0);

   helper->destroy = pfrh_destroy;
   helper->endOfRecords = pfrh_endOfRecords;
   return (PIO_HELPER*) helper;
}

void pfrh_destroy(PIO_HELPER* this)
{
   return;
}

size_t pfrh_endOfRecords(const char* buf,
			 size_t nBuf,
			 PIO_HELPER* thisBase)
{
   PIO_FIXED_RECORD_HELPER* this = (PIO_FIXED_RECORD_HELPER*) thisBase;
   size_t nRecords = nBuf / this->lrec;

   return nRecords * this->lrec;
}



/* Local Variables: */
/* tab-width: 3 */
/* End: */
