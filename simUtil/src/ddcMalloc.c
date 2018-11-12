// $Id$

#define _XOPEN_SOURCE 600
#include "ddcMalloc.h"
#include "mpiUtils.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <libgen.h>
#include <string.h>
#ifndef __APPLE__
#include <malloc.h>
#endif

#ifdef WITH_PIO
#include "pio.h"
#endif


static int addBlock(void* ptr, size_t size, char* location);
static int updateBlock(void* old_ptr, void* new_ptr, size_t size, char* location);
static int freeBlock(void* ptr);
static int findBlock(void* ptr);
static void printHeapInfo(FILE* file);
#ifdef WITH_PIO
static void printHeapInfo_pio(PFILE* file);
#endif

static int _verboseTask = -1;

typedef struct memBlock_st
{
   size_t size;
   void*  ptr;
   char   location[40];
} MEMBLOCK;

#define MAX_BLOCK 32000
static MEMBLOCK _block[MAX_BLOCK];
static int      _freeBlock[MAX_BLOCK];
static int      _nextBlock = MAX_BLOCK +1;
static size_t   _memUsed = 0;
static size_t   _peakUsed = 0;
static int      _blocksUsed = 0;
const double    b2mb=1024*1024;


void ddcMemSetVerbose(int task)
{
   _verboseTask = task;
}

void ddcMemInit(void)
{
   for (unsigned ii=0; ii<MAX_BLOCK; ++ii)
   {
      _freeBlock[ii] = ii;
      _block[ii].size = 0;
      _block[ii].ptr = NULL;
      _block[ii].location[0] = '\0';
   }
   
   _nextBlock = 0;
}

void ddcMemSummary(FILE* file)
{
   fprintf(file, 
	   "ddcMem task %d: peak=%7.2fMB current=%7.2fMB in %d blocks\n",
	   getRank(0), _peakUsed/b2mb, _memUsed/b2mb, _blocksUsed);
}

#ifdef WITH_PIO
void ddcMemSummary_pio(PFILE* file)
{
   Pprintf(file, 
	   "ddcMem task %d: peak=%7.2fMB current=%7.2fMB in %d blocks\n",
	   getRank(0), _peakUsed/b2mb, _memUsed/b2mb, _blocksUsed);
}
#endif

void ddcMemReport(FILE* file)
{
   size_t totalSize = 0;
   
   fprintf(file, "ddcMem report for task %d\n\n", getRank(0));
   fprintf(file, "Block  ptr         size          location\n");
   fprintf(file, "=======================================================================\n");
   
   for (unsigned ii=0; ii<MAX_BLOCK; ++ii)
   {
      if (_block[ii].ptr == NULL)
	 continue;
      fprintf(file, "%5d: %10p %12zuk %s\n", ii, _block[ii].ptr,
	      _block[ii].size/1024, _block[ii].location);
      totalSize += _block[ii].size;
   }
   fprintf(file, "\nTotal size = %f MB\n", totalSize/b2mb);
   fprintf(file, "Peak size = %f MB\n\n", _peakUsed/b2mb);
   printHeapInfo(file);
}

#ifdef WITH_PIO
void ddcMemReport_pio(PFILE* file)
{
   size_t totalSize = 0;
   
   Pprintf(file, "ddcMem report for task %d\n\n", getRank(0));
   Pprintf(file, "Block  ptr         size          location\n");
   Pprintf(file, "=======================================================================\n");
   
   for (unsigned ii=0; ii<MAX_BLOCK; ++ii)
   {
      if (_block[ii].ptr == NULL)
	 continue;
      Pprintf(file, "%5d: 0x%08x %12ik %s\n", ii, _block[ii].ptr,
	      _block[ii].size/1024, _block[ii].location);
      totalSize += _block[ii].size;
   }
   Pprintf(file, "\nTotal size = %f MB\n", totalSize/b2mb);
   Pprintf(file, "Peak size = %f MB\n\n", _peakUsed/b2mb);
   printHeapInfo_pio(file);
}
#endif


/** Implementation Note: Some implementations of malloc (such a purple)
 *  return a null pointer when called with a size of zero.  This is
 *  allowed by the POSIX standard.  It is entirely likely that we will
 *  call malloc with zero size since some tasks may logically have zero
 *  of some items such as fftChannels.
 *
 *  We don't want malloc to return NULL for two reasons: First, we want
 *  to use a NULL return as a sign that malloc has failed.  Second, we
 *  would like a unique pointer for each call to malloc so that we can
 *  add it to the block table.  If two entries in the block table have
 *  the same pointer they become indistinguishable when we try to free
 *  the pointer (since only the pointer is passed to ddcFree).
 *
 *  To work around this problem we check for a zero size and always
 *  allocate at least sizeof(void*).  This way a NULL return is a true
 *  error and we get a unique pointer value to add to the block table.
 */
void* _ddcMalloc(size_t size, char* location)
{
   if (size == 0)
      size = sizeof(void*);
   
   void* ptr = malloc(size);
   if (!ptr)
   {
      printf("mem: ddcMalloc failed on task %d (%zu bytes at %s)\n"
	     "                                 memUsed=%8.2f\n",
	     getRank(0), size, location, _memUsed/b2mb);
      printHeapInfo(stdout);
   }
   else
   {
      int b = addBlock(ptr, size, location);   
      if (_verboseTask == getRank(0))
      {
	 printf("mem: task %d block %d  malloc %10p  (%zu bytes at %s) total %8.2f\n",
		getRank(0), b, ptr, size, location, _memUsed/b2mb);
	 printHeapInfo(stdout);
      }
   }
   
   return ptr;
}

/** For the same reasons as explained above _ddcMalloc we check for zero
 *  size allocations and ensure at least a little memory is allocated. */
void* _ddcCalloc(size_t count, size_t size, char* location)
{
   if (count == 0)
      count = 1;
   if (size == 0)
      size = sizeof(void*);
   
   
   void* ptr = calloc(count, size);
   if (!ptr)
   {
      printf("mem: ddcCalloc failed on task %d (%zu bytes at %s\n"
	     "                                 memUsed=%8.2f\n",
	     getRank(0), size*count, location, _memUsed/b2mb);
      printHeapInfo(stdout);
   }
   else
   {
      int b = addBlock(ptr, size*count, location); 
      if (_verboseTask == getRank(0))
	 printf("mem: task %d block %d  calloc %10p  (%zu bytes at %s) total %8.2f\n",
		getRank(0), b, ptr, size*count, location, _memUsed/b2mb);
   }
   
   return ptr;
}

/** For the same reasons as explained above _ddcMalloc we check for zero
 *  size allocations and ensure at least a little memory is allocated.
 *  By POSIX, if size is zero and ptr is non-null the object is freed.*/
void* _ddcRealloc(void* ptr, size_t size, char* location)
{
   if (size == 0 && ptr == NULL)
      size = sizeof(void*);
   
   if (size == 0 && ptr != NULL)
   {
      _ddcFree(ptr, location);
      return NULL;
   }

   void* old_ptr = ptr;
   void* new_ptr = realloc(ptr, size);
   if (!new_ptr)
   {
      printf("mem: ddcRealloc failed on task %d (%zu bytes at %s)\n"
	     "                                  ptr=%10p\n"
	     "                                  memUsed=%8.2f\n",
	     getRank(0), size, location, ptr, _memUsed/b2mb);
      printHeapInfo(stdout);
   }
   else
   {
      int b = updateBlock(old_ptr, new_ptr, size, location);
      if (_verboseTask == getRank(0))
	 printf("mem: task %d block %d  realloc %10p  (%zu bytes at %s) total %8.2f\n",
		getRank(0), b, ptr, size, location, _memUsed/b2mb);
   }
   
   return new_ptr;
}

/** For the same reasons as explained above _ddcMalloc we check for zero
 *  size allocations and ensure at least a little memory is allocated. */
int _ddcMallocAligned(void** ptr, size_t alignment, size_t size, char* location)
{
   if (size == 0)
      size = sizeof(void*);
   
   int retVal = _mallocAligned(ptr, alignment, size);
   if (!*ptr)
   {
      printf("mem: ddcMallocAligned failed on task %d (%zu bytes at %s)\n"
	     "                                         memUsed=%8.2f\n",
	     getRank(0), size, location, _memUsed/b2mb);
      printHeapInfo(stdout);
   }
   else
   {
      int b = addBlock(ptr, size, location);   
      if (_verboseTask == getRank(0))
      {
	 printf("mem: task %d block %d  mallocAligned %10p  (%zu bytes at %s) total %8.2f\n",
		getRank(0), b, ptr, size, location, _memUsed/b2mb);
	 printHeapInfo(stdout);
      }
   }
   
   return retVal;
}


void  _ddcFree(void* ptr, const char* location)
{
   free(ptr);
   
   int b = freeBlock(ptr);
   
   if (_verboseTask == getRank(0))
      printf("mem:  task %d block %d  free %10p at %s total %8.2f\n",
	     getRank(0), b, ptr, location, _memUsed/b2mb);
}

char* _ddcLine(const char* file, int lineNum)
{
   static char buffer[256];
   sprintf(buffer, "%s:%d", file, lineNum);
   return buffer;
}

int addBlock(void* ptr, size_t size, char* location)
{
   assert(ptr != NULL);
   int here;
   #pragma omp critical (ddcMalloc_addBlock)
   {
      if (_nextBlock == MAX_BLOCK+1)
	 ddcMemInit();
      
      here = _freeBlock[_nextBlock];
      _block[here].ptr = ptr;
      _block[here].size = size;
      _block[here].location[0] = '\0';
      strncat(_block[here].location, basename(location),39);
      
      ++_blocksUsed;
      _memUsed += size;
      if (_memUsed > _peakUsed) _peakUsed = _memUsed;
      ++_nextBlock;
      if (_nextBlock == MAX_BLOCK)
      {
	 printf("Block storage exhausted on task %d in addBlock.\n" 
		"%s\n"
		"Try increasing MAX_BLOCK\n", getRank(0),location);
	 exit(3);
      }
   } 
   return here;
}


int updateBlock(void* old_ptr, void* new_ptr, size_t size, char* location)
{
   if (old_ptr == NULL)
      return addBlock(new_ptr, size, location);

   int here = findBlock(old_ptr);
   
   if (here < 0)
   {
      /*       printf("Error in updateBlock on task %d\n" */
      /* 	     "   old_ptr=%08x not found in block table.\n" */
      /* 	     "   new_ptr=%08x\n" */
      /* 	     "   %d bytes at %s\n", */
      /* 	     getRank(0), old_ptr, new_ptr, size, location); */
      return addBlock(new_ptr, size, location);
   }
   
   #pragma omp critical (ddcMalloc_updateBlock)
   {
      _memUsed += (size - _block[here].size);
      if (_memUsed > _peakUsed) _peakUsed = _memUsed;
      _block[here].ptr = new_ptr;
      _block[here].size = size;
      _block[here].location[0] = '\0';
      strncat(_block[here].location, basename(location),39);
   }
   return here;
}


int freeBlock(void* ptr)
{
   if (ptr == NULL)
      return -2;
   int here = findBlock(ptr);

   #pragma omp critical (ddcMalloc_freeBlock)
   {
      if (here >= 0)
      {
	 --_blocksUsed;
	 _memUsed -= _block[here].size;
	 _block[here].ptr = NULL;
	 _block[here].size = 0;
	 _block[here].location[0] = '\0';
	 --_nextBlock;
	 assert (_nextBlock >= 0);
	 _freeBlock[_nextBlock] = here;
      }
      /*    else */
      /*       printf("mem: Error on Task %d. Request to free ptr 0x%08x.\n" */
      /* 	     "     Pointer cannot be found in block list.\n", */
      /* 	     getRank(0), ptr); */
   }
   return here;
}


int findBlock(void*ptr)
{
   for (unsigned ii=0; ii<MAX_BLOCK; ++ii)
      if (_block[ii].ptr == ptr)
	 return ii;
   return -1;
}

void printHeapInfo(FILE* file)
{
#ifdef __APPLE__
 fprintf(file, "In routine printHeapInfo no mallinfo on OS X. Sorry.\n");
#else
   struct mallinfo minfo;
   minfo = mallinfo();
   fprintf(file, "mem: task %d  system heap=%fMB used=%fMB unused=%fMB\n",
	   getRank(0),
	   minfo.arena/(1024*1024.),
	   minfo.uordblks/(1024*1024.),
	   minfo.fordblks/(1024*1024.));
#endif
}

#ifdef WITH_PIO
void printHeapInfo_pio(PFILE* file)
{
#ifdef __APPLE__
   Pprintf(file, "No mallinfo on OS X.  Sorry.\n");
#else
   struct mallinfo minfo;
   minfo = mallinfo();
   Pprintf(file, "mem: task %d  system heap=%fMB used=%fMB unused=%fMB\n",
	   getRank(0),
	   minfo.arena/(1024*1024.),
	   minfo.uordblks/(1024*1024.),
	   minfo.fordblks/(1024*1024.));
#endif
}
#endif // ifdef WITH_PIO

int _mallocAligned(void** ptr, size_t alignment, size_t size)
{
#ifndef __APPLE__
   return posix_memalign(ptr, alignment, size);
#else
   *ptr = malloc(size);
   if (*ptr == NULL)
      return 1;
   return 0;
#endif
}
void freeNull(void **ptr) 
{ 
	if ( *ptr != NULL ) 
	{ free(*ptr); *ptr = NULL; }
}
