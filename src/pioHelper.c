// $Id$

#include "pioHelper.h"

#include <string.h>
#include <assert.h>

#include "ddcMalloc.h"
#include "pioFixedRecordHelper.h"
#include "pioVariableRecordHelper.h"

PIO_HELPER* pioHelperFactory(OBJECT* header)
{
   char* dataType;
   object_get(header, "datatype", &dataType, STRING, 1, "undefined");

   PIO_HELPER* helper = NULL;
   if      (strcmp(dataType, "FIXRECORDASCII") == 0)
      helper = pfrh_create(header);
   else if (strcmp(dataType, "FIXRECORDBINARY") == 0)
      helper = pfrh_create(header);
   else if (strcmp(dataType, "VARRECORDASCII") == 0)
      helper = pvrah_create(header);
   else if (strcmp(dataType, "VARRECORDBINARY") == 0)
      helper = pvrbh_create(header);

   // catch attempt to read unknown datatype.
   assert(helper != NULL);
   
   ddcFree(dataType);
   return helper;
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
