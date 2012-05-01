#ifndef FASTBARRIER_NOSYNC__
#define FASTBARRIER_NOSYNC__

#ifdef FAST_BARRIER_HH

#define PER_SQUAD_BARRIER

#ifdef FAST_BARRIER_BGQ_HH
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
  }
  // prevent speculation
  __isync();
}

#else
  /* No fast bgq barrier, use regular synchronizing fast barrier 
     for the per squad barriers.
  */

__INLINE__ void L2_Barrier_nosync_Arrive(
  L2_Barrier_t *b,        /* global barrier */
  L2_BarrierHandle_t *h,  /* barrier handle private to this thread */
  int eventNum)           /* number of arrival events */
{
  L2_BarrierWithSync_Arrive(b,h,eventNum);
}
__INLINE__ void L2_Barrier_nosync_WaitAndReset(
  L2_Barrier_t *b,       /* global barrier */
  L2_BarrierHandle_t *h, /* barrier handle private to this thread */
  int eventNum)          /* number of arrival events */
{
  L2_BarrierWithSync_WaitAndReset(*b,h,eventNum);
}


#endif /* Fast bgq barrier if-statement */

#endif /* Fast barrier if-statement */
#endif /* Multiple inclusion protection */
