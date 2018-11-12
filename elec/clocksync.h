#ifndef CLOCKSYNC__
#define CLOCKSYNC__

#include "cs_gettime.h"

#ifdef __cplusplus
extern "C" {
#endif

  /*
    Local time is defined by whatever cs_gettime() returns.
    cs_gettime() returns a double precision number, and is
    defined in cs_gettime.c. All times are measured in the
    same unit as returned by cs_gettime(), whever that may
    be.

    This subroutine computes an offset tc to local time, such
    that cs_gettime()+tc is a globally synchronized time. The
    task with mpi rank 0 has tc = 0.0 by definition.

    Input arguments:
    np -- number of MPI ranks
    pid -- MPI rank of current task.
    nmsg -- Number of messages to send/receive to calculate
            a time correction. Accuracy of synchronization
	    is proportional to 1/sqrt(nmsg).
    gettime_cost:
      If gettime_cost is NULL or gettime_cost[0] < 0.0, then
        the time to make a call to cs_gettime() is measured.
      otherwise,
        the time to make a call to cs_gettime() is assumed
	to be gettime_cost[0].
    corr_std:
      If not NULL, corr_std[0] will contain the standard
      deviation of the measured time correction on return.
   */

  double cs_clocksync(int np,int pid,int nmsg,
		      double gettime_cost[1],
		      double corr_std[1]);

#ifdef __cplusplus
}
#endif


#endif
