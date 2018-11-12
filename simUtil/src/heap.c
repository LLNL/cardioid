/* $Id$ */

#include "heap.h"
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <libgen.h>
#include "mpiUtils.h"
#include "ddcMalloc.h"
#include "object.h"

#define MAX_HEAP_BLOCKS 64
extern FILE *ddcfile; 

// Is there any reason we shouldn't quad align all of the heap
// allocations on every platform?  I can't think of any.
#define QUAD_WORD_ALIGN_HEAP
static const int alignOffset=16;

typedef struct heapBlock_st
{
	char* blockStart;
	size_t blockSize;
	char location[40];
} HEAP_BLOCK;


static void* heapBuffer=NULL;
static size_t _heapCapacity = 0;
static size_t _heapClaimed = 0;
static int _nBlocks=0;
static HEAP_BLOCK _blocks[MAX_HEAP_BLOCKS];
static int _openBlock = -1;
static size_t _heapHighWaterMark = 0;

static void blockOpenError(const char* file, int line);
static void tooManyBlocksError(const char* file, int line);
static void heapTooSmallError(size_t request, const char* file, int line);
static void heapFreeOpenBlockError(const char* file, int line);
static void heapFreeOrderError(const char* file, int line);
static void heapAlreadyAllocatedWarning(void);
static void heapAlreadyClosedError(const char* file, int line);


HEAP* heap_init(void* parent, char* name)
{
	HEAP* heapObj = (HEAP*)object_initialize(name, "HEAP", sizeof(HEAP));
	heapObj->parent=parent; 
	int size;
	object_get((OBJECT *) heapObj, "size", &size, INT, 1, "1");

	/* The size in bytes can be large, so make sure to compute as size_t */ {
	  size_t l_size = size;
	  l_size *= (size_t) (1024*1024);
	  heap_allocate(l_size);
	}
	_nBlocks = 0;
	_openBlock = -1;
	for (unsigned ii=0; ii<MAX_HEAP_BLOCKS; ++ii)
	{
		_blocks[ii].blockStart = NULL;
		_blocks[ii].blockSize = 0;
		_blocks[ii].location[0] = '\0';
	}
	
	return heapObj; 
}
void heap_start(int size)
{

  /* The size in bytes can be large, so make sure to compute as size_t */ {
	  size_t l_size = size;
	  l_size *= (size_t) (1024*1024);
	  heap_allocate(l_size);
	}
	_nBlocks = 0;
	_openBlock = -1;
	for (unsigned ii=0; ii<MAX_HEAP_BLOCKS; ++ii)
	{
		_blocks[ii].blockStart = NULL;
		_blocks[ii].blockSize = 0;
		_blocks[ii].location[0] = '\0';
	}
}

void* _heapGet(unsigned* blockNumber, char* file, int line)
{
	if (_openBlock != -1)
		blockOpenError(file, line);

	_openBlock = _nBlocks;
	*blockNumber= _nBlocks;
	//fprintf(ddcfile,"ddt %d:  Opening block %d (%s:%d)\n", getRank(-1), _nBlocks, file, line); fflush(ddcfile);
	
	void* blockStart = ((char*)heapBuffer)+_heapClaimed;
	_blocks[_nBlocks].blockStart = blockStart;
	_blocks[_nBlocks].blockSize  = _heapCapacity-_heapClaimed;
	snprintf(_blocks[_nBlocks].location, 39, "%s:%d", basename(file), line);
	++_nBlocks;
	if (_nBlocks >= MAX_HEAP_BLOCKS)
		tooManyBlocksError(file, line);


	
	return blockStart;
}
void slave_heap_init(char *data) { heap_init(NULL,data); }
void _heapTestAvailable(unsigned blockNumber,  size_t size, const char* file, int line)
{
	size_t heapNeeded = _heapClaimed + size; 
	if (heapNeeded > _heapCapacity) heapTooSmallError(size, file, line);
}
size_t heapAvailable(void) {
  return _heapCapacity - _heapClaimed - 64 /* In case 32 bytes are lost on either side for alignment reasons */;
}
void _heapEndBlock(unsigned blockNumber, size_t size, const char* file, int line)
{
	//fprintf(ddcfile,"ddt %d:  Ending block %d (%s:%d), size=%d\n", getRank(-1), blockNumber, file, line, size); fflush(ddcfile); 
	if (_openBlock != (int)blockNumber) heapAlreadyClosedError(file, line);
	
	_openBlock = -1;

#ifdef QUAD_WORD_ALIGN_HEAP
	size += (alignOffset - (size%alignOffset) ) % alignOffset;
#endif


	assert((int)blockNumber == _nBlocks-1);
	
	_heapClaimed += size; 
	_blocks[_nBlocks-1].blockSize = size; 
	if (_heapHighWaterMark < _heapClaimed) _heapHighWaterMark = _heapClaimed; 
	if (_heapHighWaterMark > _heapCapacity) heapTooSmallError(size, file, line);
}
void _heapFree(unsigned blockNumber, const char* file, int line)
{
 	//fprintf(ddcfile,"ddt %d:  Free block %d (%s:%d) size=%d, claimed=%d\n", getRank(-1), blockNumber, file, line, _blocks[_nBlocks-1].blockSize, _heapClaimed); fflush(ddcfile); 
	if (_openBlock != -1)
		heapFreeOpenBlockError(file, line);
	if ((int)blockNumber != _nBlocks-1)
		heapFreeOrderError(file, line);

	if (heapBuffer == NULL || _nBlocks == 0)
		return ;

	--_nBlocks;
	assert(_nBlocks == (int)blockNumber);
	_heapClaimed -= _blocks[_nBlocks].blockSize;

	_blocks[_nBlocks].blockStart = NULL;
	_blocks[_nBlocks].blockSize = 0;
	_blocks[_nBlocks].location[0] = '\0';
}

void heap_allocate(size_t n)
{
	if (heapBuffer != NULL)
	{
		heapAlreadyAllocatedWarning();
		return; 
	}
	
	_heapHighWaterMark = 0; 
	_heapClaimed = 0;
	_heapCapacity = n; 
#ifdef  QUAD_WORD_ALIGN_HEAP
	_heapCapacity +=  (alignOffset - (_heapCapacity % alignOffset) ) % alignOffset;
	ddcMallocAligned((void **)&heapBuffer, alignOffset, _heapCapacity);
#else 
	heapBuffer = ddcMalloc(n); 
#endif 
	_nBlocks=0; 
	
}
size_t heapSize(void)
{
	return _heapCapacity;
}

size_t heapBlockSize(unsigned block)
{
   if ((int)block >= _nBlocks)
      return 0;
	return _blocks[block].blockSize;
}

void heapSummary(FILE *file)
{
	fprintf(file, "Scratch heap capacity: %f Mbytes,  high water mark: %f Mbytes \n",
			  _heapCapacity/(1024*1024.0), _heapHighWaterMark/(1024*1024.0));
}

#ifdef WITH_PIO
void heapSummary_pio(PFILE *file)
{
	Pprintf(file, "Scratch heap capacity: %f Mbytes,  high water mark: %f Mbytes \n",
			  _heapCapacity/(1024*1024.0), _heapHighWaterMark/(1024*1024.0));
}
#endif

void heap_deallocate(void)
{
	ddcFree(heapBuffer); 
}


void blockOpenError(const char* file, int line)
{
	printf("Scratch Heap ERROR:\n"
			 "   heapGet called with block already open.\n"
			 "   current call from %s:%d\n"
			 "   open block created from %s\n",
			 file, line, _blocks[_nBlocks-1].location);
	abortAll(-1);
}

void tooManyBlocksError(const char* file, int line)
{
	printf("Scratch Heap ERROR:\n"
			 "   Too many blocks.  _nBlocks = %d\n"
			 "   current call from %s:%d\n",
			 _nBlocks, file, line);
	abortAll(-1);
}

void heapTooSmallError(size_t request, const char* file, int line)
{
	printf("Scratch Heap ERROR:\n"
			 "   heap too small\n"
			 "   current call from %s:%d\n"
			 "   _heapCapacity = %f Mbytes\n"
			 "   _heapClaimed = %f Mbytes\n"
			 "   request = %f Mbytes\n",
			 file, line, _heapCapacity/(1024*1024.0),
			 _heapClaimed/(1024*1024.0),
			 request/(1024*1024.0));
	abortAll(-1);
}

void heapFreeOpenBlockError(const char* file, int line)
{
	printf("Scratch Heap ERROR:\n"
			 "   Attempt to free open block from %s:%d\n",
			 file, line);
	abortAll(-1);
}

void heapFreeOrderError(const char* file, int line)
{
	printf("Scratch Heap ERROR:\n"
			 "   Attempt to free blocks out of order %s:%d\n",
			 file, line);
	abortAll(-1);
}

void heapAlreadyAllocatedWarning(void)
{
}

void heapAlreadyClosedError(const char* file, int line)
{
	printf("Scratch Heap ERROR:\n"
			 "   Attempt to call heapEndBlock on a block that is already closed\n"
			 "   %s:%d\n",
			 file, line);
	abortAll(-1);
}




/* Local Variables: */
/* tab-width: 3 */
/* End: */
