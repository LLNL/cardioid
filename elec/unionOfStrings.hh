#ifndef UNION_OF_STRINGS_HH
#define UNION_OF_STRINGS_HH

#include <mpi.h>
#include <vector>
#include <string>

void unionOfStrings(std::vector<std::string>& words, MPI_Comm comm, int tag);

#endif
