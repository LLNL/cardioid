/* $Id$ */
#ifndef PIO_H
#define PIO_H

#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE 1
#define _LARGE_FILE
#include <stdio.h>
#include <mpi.h>
#include "object.h"
#include "pioHelper.h"

/**
 *  LIMITATIONS:
 *
 *  The biggest limitation of pio is that it does not support a
 *  streaming I/O model.  The entire contents of the pfile must be
 *  present in memory at the start/end of each read/write.  The data is
 *  distributed across the tasks however.  Most tasks will create a
 *  buffer large enough to store all of the data that will be read or
 *  written by that task.  The reader/writer tasks will need to store
 *  their own data plus a buffer big enough to send/recv the data for
 *  one additional task.  For memory bound applications with large I/O
 *  requirements there may be insufficient memory to create the
 *  necessary buffers.
 */

/** To be readable by pio a file must have and OBJECT style file header
 *	 with the following properties:
 *
 *  #  The header must be followed by two consecutive newline to
 *     delimit the end of the header.  The header cannot contain
 *     two consecutive newlines.
 *  #  The name and class of the OBJECT are arbitrary.
 *  #  The header must provide a value for the nfiles keyword
 *     indicating the number of physical files that are part of
 *     the pfile.
 *  #  The header must provide a value for the datatype keyword
 *     to indicate how the data formatted.  This information is
 *     used to create a PIO_HELPER object that allows pio to
 *     distribute data correctly to tasks in an I/O group.
 */

#ifdef __cplusplus
extern "C" {
#endif

/* if you change PIO_ENUMS please update the definition of PioNames in
 * pio.c to keep it in sync. */
enum PIO_ENUMS { PIO_NONE, SPLIT, FIXRECORDASCII, FIXRECORDBINARY,
					  VARRECORDASCII, VARRECORDBINARY, CRC32};
extern char* PioNames[];

typedef unsigned long long pio_long64;
	
typedef struct pfile_st
{
	int id;
	int size;
	int ngroup;
	int nfiles;
	int sizegroup;
	int io_id;
	int groupToHandle;
	int beanCounter;
	int doneTag;
	int goTag;
	int msgTag;
	unsigned datatype;
	unsigned recordLength;
	unsigned headerLength;
	unsigned nfields;
	int recordIndex; // used only by check_line.c
	pio_long64 numberRecords; 
	pio_long64 filesize; 
	int checksum; 
	FILE* file;
	FILE** readFile; /* when reading we may have more than 1 file per task. */
	pio_long64* nBytesInFile; /* number of bytes in each read file */
	char *buf, *name, *mode, *field_names,*field_types,*field_units,*misc_info;
	unsigned pio_buf_blk;
   size_t bufsize, bufcapacity, bufpos;
	OBJECT* headerObject;
	PIO_HELPER* helper;
	MPI_Comm comm;
} PFILE;

/** When writing data pio requires a scratch heap large enough to hold
 *  the output message plus the mpi buffer to recive that message on the
 *  writer task.  So, allocate a scratch heap at least twice the size of
 *  the maximum data written per task.
 *
 *  When reading data pio does not require a scratch heap (it gets
 *  memory from malloc).
 */

PFILE* Popen(const char *filename, const char *mode, MPI_Comm comm);
int Pclose(PFILE* file);
char* Pfgets(char* string, int size, PFILE* file);
void PioSet(PFILE*file, const char *string, ... );
void Pget(PFILE*file, const char *string, void *ptr);
size_t Pwrite(const void *ptr, size_t size, size_t nmemb, PFILE*file);
size_t Pread(void* ptr, size_t size, size_t nitems, PFILE* file);
void slave_Pio_setNumWriteFiles(int* nWriteFiles);
void Pio_setNumWriteFiles(int nWriteFiles);
void PioReserve(PFILE* file, size_t size);
// What should be the behavior of Pprintf when the line is too long for
// the 1024 character buffer that is allocated?  Right now the line is
// silently truncated.
int Pprintf(PFILE*file, const char *fmt, ...);

#ifdef __cplusplus
}
#endif

#endif

/* Local Variables: */
/* tab-width: 3 */
/* End: */
