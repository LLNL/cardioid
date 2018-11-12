#ifndef READ_PIO_FILE_HH
#define READ_PIO_FILE_HH

struct pfile_st;
class BucketOfBits;

/** Caller must delete returned pointer */
BucketOfBits* readPioFile(pfile_st* file);

#endif
