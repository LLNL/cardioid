#ifndef ANATOMY_READER_H
#define ANATOMY_READER_H

#include <string>
#include <set>
#include <mpi.h>

#include "AnatomyCell.hh"

struct pfile_st;
class Anatomy;
class BucketOfBits;

/** Caller must delete returned pointer */
BucketOfBits* readAnatomy(const std::string& filename, MPI_Comm comm,
                          Anatomy& anatomy, const std::set<int>& typeSet);

#endif
