#include "intQueue.h"
#include <assert.h>

#include "ddcMalloc.h"

#define MAX(A,B) ((A) > (B) ? (A) : (B))

INT_QUEUE* intQueue_init(unsigned capacity)
{
   INT_QUEUE* this = ddcMalloc(sizeof(INT_QUEUE));
   this->begin = 0;
   this->end = 0;
   this->size = 0;
   this->capacity = capacity;
   this->data = ddcMalloc(capacity*sizeof(int));

   return this;
}

void intQueue_destroy(INT_QUEUE* this)
{
   ddcFree(this->data);
}

void intQueue_push(INT_QUEUE* this, int item)
{
   if (this->size == this->capacity)
      intQueue_reserve(this, 2*this->capacity);

   this->data[this->end] = item;

   ++this->end;
   if (this->end == this->capacity)
      this->end = 0;
   ++this->size;
}

int intQueue_pop(INT_QUEUE* this)
{
   assert(this->size > 0);
   int tmp = this->begin;
   ++this->begin;
   --this->size;
   if (this->begin == this->capacity)
      this->begin =0;
   
   return this->data[tmp];
}

unsigned intQueue_size(const INT_QUEUE* const this)
{
   return this->size;
}

unsigned intQueue_capacity(const INT_QUEUE* const this)
{
   return this->capacity;
}


/** Ensure the queue has the capcity to store at least capacity items.
 * In the case of a wrapped queue, (end < begin) move the data so that
 * end > begin.*/
void intQueue_reserve(INT_QUEUE* this, unsigned capacity)
{
   if (capacity <= this->capacity)
      return;

   if (this->size == 0 || this->end > this->begin)
   {
      this->data = ddcRealloc(this->data, capacity*sizeof(int));
   }
   else
   {
      capacity = MAX(this->capacity+this->end, capacity);
      this->data = ddcRealloc(this->data, capacity*sizeof(int));
      for (unsigned ii=0; ii<this->end; ++ii)
	 this->data[this->capacity+ii] = this->data[ii];
      this->end += this->capacity;
   }
   this->capacity = capacity;
   return;   
}

void intQueue_clear(INT_QUEUE* this)
{
   this->begin=0;
   this->end=0;
   this->size=0;
}
