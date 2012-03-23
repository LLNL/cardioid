#ifndef CHECKPOINT_IO_HH
#define CHECKPOINT_IO_HH

#include <mpi.h>
#include <string>
class Simulate;

void writeCheckpoint(const Simulate& sim, MPI_Comm comm);
void readCheckpoint(const std::string& filename, Simulate& sim, MPI_Comm comm);
#endif
