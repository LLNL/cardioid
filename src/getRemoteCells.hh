#ifndef GET_REMOTE_CELLS_HH
#define GET_REMOTE_CELLS_HH

#include <mpi.h>
#include <string>

class Simulate;

void getRemoteCells(Simulate& sim, const std::string& name, MPI_Comm comm);

#endif
