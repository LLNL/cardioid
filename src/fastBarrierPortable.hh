#ifndef FAST_BARRIER_PORTABLE_HH
#define FAST_BARRIER_PORTABLE_HH

#include <stdint.h>
#include <cstdlib>

#ifndef FAST_BARRIER_HH
#error "Do not #include fastBarrierPortable.hh.  #include fastBarrier.hh instead"
#endif


struct  L2_Barrier_t
{
   volatile  uint64_t start;  /*!< Thread count at start of current round. */
   volatile  uint64_t count;  /*!< Current thread count. */
};


struct L2_BarrierHandle_t
{
  uint64_t localStart; // local (private start)
  volatile uint64_t *localCountPtr; // encoded fetch and inc address
};


inline void L2_BarrierWithSync_Init(L2_Barrier_t* b)
{
   b->start = b->count = 0;
} 

inline void L2_BarrierWithSync_InitInThread(
   L2_Barrier_t *b,       /* global barrier */
   L2_BarrierHandle_t *h) /* barrier handle private to this thread */
{
   h->localStart = 0;
   h->localCountPtr = &b->count;
}


inline void L2_BarrierWithSync_Arrive(
   L2_Barrier_t *b,        /* global barrier */
   L2_BarrierHandle_t *h,  /* barrier handle private to this thread */
   int eventNum)           /* number of arrival events */
{
   #pragma omp atomic
   (*(h->localCountPtr))++;
//   uint64_t current = count;
   uint64_t current = *h->localCountPtr;
   // if reached target, update barrier's start
   uint64_t target = h->localStart + eventNum;
   if (current == target)
   {
      b->start = current;  // advance to next round
   }
}

inline void L2_BarrierWithSync_WaitAndReset(
   L2_Barrier_t *b,       /* global barrier */
   L2_BarrierHandle_t *h, /* barrier handle private to this thread */
   int eventNum)          /* number of arrival events */
{
   // compute target from local start
   uint64_t target = h->localStart + eventNum;
   // advance local start
   h->localStart = target;
   // wait until barrier's start is advanced
   while (b->start < target) ;
}

inline void L2_BarrierWithSync_Reset(
   L2_Barrier_t *b,       /* global barrier */
   L2_BarrierHandle_t *h, /* barrier handle private to this thread */
   int eventNum)          /* number of arrival events */
{
   // compute target from local start
   uint64_t target = h->localStart + eventNum;
   // advance local start
   h->localStart = target;
}

inline void L2_BarrierWithSync_Barrier(
   L2_Barrier_t *b,       /* global barrier */
   L2_BarrierHandle_t *h, /* barrier handle private to this thread */
   int eventNum)          /* number of arrival events */
{
   L2_BarrierWithSync_Arrive(b, h, eventNum);
   L2_BarrierWithSync_WaitAndReset(b, h, eventNum);
}

/** Call this before the parallel region.  Replaces call to Init.
 * Caller is responsible to free the returned pointer. */
inline L2_Barrier_t* L2_BarrierWithSync_InitShared()
{
   L2_Barrier_t* bb = (L2_Barrier_t*) malloc(sizeof(L2_Barrier_t));
   L2_BarrierWithSync_Init(bb);
   return bb;
}

   


#endif
