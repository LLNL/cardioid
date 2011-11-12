// $Id:

#include "lessThan.h"

int unsignedLess(const void* a, const void* b)
{
   unsigned u1 = *((unsigned*) a);
   unsigned u2 = *((unsigned*) b);
   if (u1<u2) return -1;
   if (u1==u2) return 0;
   return 1;
}

int intLess(const void* a, const void* b)
{
   int u1 = *((int*) a);
   int u2 = *((int*) b);
   if (u1<u2) return -1;
   if (u1==u2) return 0;
   return 1;
}



/* Local Variables: */
/* tab-width: 3 */
/* End: */
