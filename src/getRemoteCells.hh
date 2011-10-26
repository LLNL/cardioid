#ifndef GET_REMOTE_CELLS_HH
#define GET_REMOTE_CELLS_HH

#include <mpi.h>

class Simulate;

void getRemoteCells(Simulate& sim, MPI_Comm comm);

#endif
