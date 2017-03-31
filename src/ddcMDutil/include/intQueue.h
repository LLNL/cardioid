#ifndef INT_QUEUE_H
#define INT_QUEUE_H

#ifdef __cplusplus
extern "C" 
{
#endif

typedef struct IntQueue_st
{
   unsigned begin;
   unsigned end;
   unsigned size;
   unsigned capacity;
   int* data;
}
INT_QUEUE;

/** Implements a first-in first-out queue of integers with push/pop
 *  interface.  The needed memory is automatically allocated.  */

INT_QUEUE* intQueue_init(unsigned capacity);
void intQueue_destroy(INT_QUEUE* intQueue);

void     intQueue_push(INT_QUEUE* intQueue, int item);
int      intQueue_pop(INT_QUEUE* intQueue);
unsigned intQueue_size(const INT_QUEUE* const intQueue);
unsigned intQueue_capacity(const INT_QUEUE* const intQueue);
void     intQueue_reserve(INT_QUEUE* intQueue, unsigned capacity);
void     intQueue_clear(INT_QUEUE* intQueue);

#ifdef __cplusplus
}
#endif

#endif
