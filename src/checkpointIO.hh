#ifndef CHECKPOINT_IO_HH
#define CHECKPOINT_IO_HH

#include <mpi.h>
class Simulate;

void writeCheckpoint(const Simulate& sim, MPI_Comm comm);
void readCheckpoint(Simulate& sim, MPI_Comm comm);
#endif
