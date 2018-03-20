// $Id$

#ifndef DDCMALLOC_H
#define DDCMALLOC_H

#include <stdlib.h>
#include <stdio.h>

struct pfile_st; // avoid including pio.h

#ifdef __cplusplus
extern "C"{
#endif


#ifdef FAST_MALLOC
// Even in the FAST_MALLOC case ddcMallocAlign does not call
// posix_memalign directly.  This is because we do not want to define
// _XOPEN_SOURCE in this header file and drag it into practically the
// entire code base.
#define ddcMalloc(size)         malloc(size)
#define ddcCalloc(count, size)  calloc(count, size)
#define ddcRealloc(ptr, size)   realloc(ptr, size)
#define ddcMallocAligned(ptr, alignment, size) \
   _mallocAligned(ptr, alignment, size)
#define ddcFree(ptr)            free(ptr)

#else

#define ddcMalloc(size) _ddcMalloc(size, _ddcLine(__FILE__, __LINE__))
#define ddcCalloc(count, size)  \
   _ddcCalloc(count, size, _ddcLine(__FILE__, __LINE__))
#define ddcRealloc(ptr, size) \
   _ddcRealloc(ptr, size, _ddcLine(__FILE__, __LINE__))
#define ddcMallocAligned(ptr, alignment, size) \
   _ddcMallocAligned(ptr, alignment, size, _ddcLine(__FILE__, __LINE__))
#define ddcFree(ptr) _ddcFree(ptr, _ddcLine(__FILE__, __LINE__))

#endif

void ddcMemSetVerbose(int task);
void ddcMemInit(void);
void ddcMemSummary(FILE*);
void ddcMemReport(FILE*);

void ddcMemSummary_pio(struct pfile_st*);
void ddcMemReport_pio(struct pfile_st*);
void freeNull(void **ptr) ;
   
void* _ddcMalloc(size_t size, char* location);
void* _ddcCalloc(size_t count, size_t size, char* location);
void* _ddcRealloc(void* ptr, size_t size, char* location);
int   _mallocAligned   (void** ptr, size_t alignment, size_t size);
int   _ddcMallocAligned(void** ptr, size_t alignment, size_t size, char* location);  
void  _ddcFree(void* ptr, const char* location);

char* _ddcLine(const char* file, int lineNum);

#ifdef __cplusplus
}
#endif

#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
