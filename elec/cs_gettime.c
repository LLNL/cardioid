#include "cs_gettime.h"

#ifdef BGQ

#include <bgpm/include/bgpm.h>

static unsigned long long int tfirst = 0;

double cs_gettime(void) {
  return (double) (long long int) (GetTimeBase() - tfirst);
}

#else

#include <time.h>
#include <sys/time.h>

static time_t tfirst = 0;

double cs_gettime(void) {
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return (tv.tv_sec - tfirst) + 1e-6*tv.tv_usec;
}

#endif

void cs_inittime(void) {
  tfirst = cs_gettime();
}

double cs_getlocaloffset(void) {
  return tfirst;
}
