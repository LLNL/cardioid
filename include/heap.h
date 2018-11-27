/* $Id$  */
#ifndef HEAP_H
#define HEAP_H

#include <stdio.h>
#ifdef WITH_PIO
#include "pio.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** Implements a "scratch heap" memory pool.  The idea behind this pool
 *  is to provide a relatively large area of memory that is available
 *  for transient calculations that persist only within a given scope.
 *  Such memory could be obtained from malloc, but this can lead to
 *  memory fragmentation especially when the amount of memory needed
 *  isn't known in advance and realloc must be called to grow an array
 *  initially allocated too small.
 *
 *  A pointer to the heap is always in one of two states: open or
 *  closed.  The heapGet function returns a pointer to an open block.
 *  An open block is allowed to access memory through the end of the
 *  scratch heap.  The heapEndBlock function converts an open block to a
 *  closed block.  Closed blocks are expected to access only memory
 *  within the size limit passed to heapEndBlock.
 *
 *  heapGet returns a pointer either to the start of the heap (if there
 *  are no closed block) or to the next valid address after the end of
 *  the last closed block.  Modules can request multiple consecutive
 *  blocks of memory within the scratch heap by repeatedly calling
 *  heapGet followed by heapEndBlock.
 *
 *  Clearly all sorts of problems will occur if any module attempts to
 *  use more memory than the scratch heap provides or if two modules
 *  attempt to use the same memory areas within the sratch heap.  To
 *  help detect such problems this implementation will enforce the
 *  following requirements:
 *
 *  1.  Only a single block can be "open" at any time.
 *  2.  heapEndBlock cannot request more space than is available in
 *      scratch heap (excluding any existing closed blocks).
 *  3.  heapFree can be called only on "closed" blocks.
 *  4.  Every call to heapGet must have a matching call to heapFree.
 *  5.  Calls to heapFree  must occur in reverse sequence of the
 *      heapGet calls.  I.e, heapFree can only be called for
 *      closed block at the end of the scratch heap. 
 *
 *
 *  Taken together, these three requirements impose the following
 *  required conventions for callers:
 *
 *  1.  heapEndBlock *must* be called for every block.  The call may be
 *      deferred until the caller knows how much space is needed, but
 *      heapEndBlock must be called even if it is immediately before
 *      heapFree.
 *  2.  Callers must call heapFree on each block they request.
 *  3.  Calls to heapGet and heapFree must nest properly at all
 *      levels of the application.  This is intended to enforce the
 *      idea that the scratch heap is temporary storage that must be
 *      released when a module goes out of scope.  It will also prevent
 *      fragmentation of the heap.
 */
   
typedef struct heap_st
{
   char *name;      /* name of the system */
   char *objclass;
   char *value;
   char *type;      /* type label */
   void  *parent;
   
} HEAP;

HEAP *heap_init(void* parent, char* name);
void heap_start(int size);
void  heap_allocate(size_t size);
void  heap_deallocate(void);
void slave_heap_init(char *data);

/** Pass a non-zero flag to active verbose tracing of heap operations */
void heap_setVerbse(int flag);
   
#define heapGet(blockNumber)  _heapGet(blockNumber, __FILE__, __LINE__)
#define heapFree(blockNumber) _heapFree(blockNumber, __FILE__, __LINE__)
#define heapEndBlock(blockNumber, blockSize) _heapEndBlock(blockNumber, blockSize, __FILE__, __LINE__)
#define heapTestAvailable(blockNumber, blockSize) _heapTestAvailable(blockNumber, blockSize, __FILE__, __LINE__)
   
void* _heapGet(unsigned *blockNumber, char* file, int line);
void  _heapFree(unsigned blockNumber, const char* file, int line);
void  _heapEndBlock(unsigned blockNumber, size_t size, const char* file, int line);
void  _heapTestAvailable(unsigned blockNumber, size_t size, const char* file, int line);
size_t heapAvailable(void);

/** Returns total size of heap. (I.e., the size that was allocated.) */
size_t  heapSize(void);

/** Returns the size of the specified block.  If heapEndBlock has not
 * been called for the specified block then the size is computed to the
 * end of the heap.  Returns 0 if the specified block has not been
 * created by heapGet (or if heapFree has been called on the block).*/
size_t heapBlockSize(unsigned blockNumber);

void heapReport(FILE* file);
void heapSummary(FILE* file);

#ifdef WITH_PIO
   void heapReport_pio(PFILE* file);
   void heapSummary_pio(PFILE* file);
#endif

#ifdef __cplusplus
}
#endif

#endif

// Questions / ToDo:
// 1.  Should the verbose keyword instead take a string that would
//     be interprested as a regular expression that would match
//     task numbers?  How can we use task numbers without introducing
//     MPI into this module?
// 2.  How concerned should we be about the fact that basename isn't
//     re-entrant?
// 3.  Code for verbose tracking isn't written
// 4.  Error messages don't report task number.
// 5.  Better way to die than MPI_Abort
// 6.  heapAlreadyAllocatedWarning doesn't do anything.




/* Local Variables: */
/* tab-width: 3 */
/* End: */
