/* $Id$ */ 
#ifndef IOUTILS_H
#define IOUTILS_H

#include <sys/stat.h>  // brings in flags that are passed as type arg to
		       // filetest function.

#ifndef S_IFDIR
#define S_IFDIR	__S_IFDIR
#endif

#ifndef S_IFREG
#define S_IFREG	__S_IFREG
#endif

#ifdef __cplusplus
extern "C"{
#endif

enum IO_FIELDS {IOF_CHECKSUM, IOF_ID, IOF_PINFO,
		IOF_RX, IOF_RY, IOF_RZ,
		IOF_VX, IOF_VY, IOF_VZ,
		IOF_NONE};

int filetest(const char* filename, int type);
int DirTestCreate(const char *dirname);


void endianSwap(void* data, int size);

void ioUtils_setSwap(unsigned endianKey);

/** Be sure to call ioUtils_setSwap before calling mkInt or mkDouble. */
unsigned long long mkInt(const unsigned char* data, const char* fieldType);
double mkDouble(const unsigned char* data, const char* fieldType);
double mkValue(const unsigned char* data, const char* fieldType);

/** Replacement for memcpy that makes no attempt to be clever.  On BG/L
 * using memcpy to copy floating point variables to non-aligned
 * addresses causes exceptions to be thrown.  Logging these exceptions
 * to the RAS database at the end of a run can take over 1 hour on large
 * partitions. */
void copyBytes(void* out, const void* in, unsigned size);

unsigned long long
bFieldUnpack(const unsigned char* buf, const char* fieldType);
void bFieldPack(unsigned char* buf, unsigned size, unsigned long long in);
unsigned bFieldSize(unsigned long long i);

unsigned* makeOffsets(char** fieldType, unsigned nFields);
enum IO_FIELDS* makeFieldMap(char** fieldNames, unsigned nFields);
unsigned parseFieldSizes(unsigned* offset, const char* field_types);
/** Returns the number of requestedNames or -1 if a requested name
 *  cannot be found among the inputNames.
 *
 *  The caller must allocate the the offset and type arrays.
 */
int makeOffsetsAndTypes(unsigned nInput, char** inputNames,
			char** inputTypes, const char* requestedNames,
			unsigned* offset, char** type);
#ifdef __cplusplus
}
#endif


#endif // #ifndef IOUTILS_H


/* Local Variables: */
/* tab-width: 3 */
/* End: */
