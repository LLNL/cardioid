#ifndef FASTBARRIER_NOSYNC_HH
#define FASTBARRIER_NOSYNC_HH


#ifdef BGQ
#define PER_SQUAD_BARRIER
  /* This barrier code is correct for barriers executing
     entirely within one core.
  */

__INLINE__ void L2_Barrier_nosync_Arrive(
  L2_Barrier_t *b,        /* global barrier */
  L2_BarrierHandle_t *h,  /* barrier handle private to this thread */
  int eventNum)           /* number of arrival events */
{
  // memory sync
  //__lwsync();
  // fetch and increment
  uint64_t count =  *(h->localCountPtr);
  uint64_t current = count + 1;
  // if reached target, update barrier's start
  uint64_t target = h->localStart + eventNum;
  if (current == target) {
      b->start = current;  // advance to next round
  }
}
__INLINE__ void L2_Barrier_nosync_WaitAndReset(
  L2_Barrier_t *b,       /* global barrier */
  L2_BarrierHandle_t *h, /* barrier handle private to this thread */
  int eventNum)          /* number of arrival events */
{
  // compute target from local start
  uint64_t target = h->localStart + eventNum;
  // advance local start
  h->localStart = target;
  // wait until barrier's start is advanced
  if (b->start < target) {
    // must wait

    if(0) { /* Sleep loop , rumorously causes SPI error */
      while (1) {
	// load start and reserve
	uint64_t remoteVal = __ldarx((volatile long *) &b->start);
	if (remoteVal >= target) {
	  // just witnessed the value to be fine
	  break;
	  // wait with sleep (may be awoken spuriously)
	  ppc_waitrsv();
	}
      }
    } else { /* Spin loop */
      while(b->start < target) {}
    }
  }
  // prevent speculation
  __isync();
}
#endif  // #ifdef BGQ
#endif  // #ifndef FASTBARRIER_NOSYNC_HH

