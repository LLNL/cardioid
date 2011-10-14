#ifndef INITIALIZE_ANATOMY_HH
#define INITIALIZE_ANATOMY_HH

#include <mpi.h>
#include <string>

class Simulate;

void initializeAnatomy(Simulate& sim, const std::string& name, MPI_Comm comm);

#endif
