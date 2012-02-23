#ifndef WRITE_CHECKPOINT_HH
#define WRITE_CHECKPOINT_HH

#include <string>
#include <mpi.h>
class Simulate;

void writeCheckpoint(const Simulate& sim, const std::string& filename, MPI_Comm comm);

#endif
