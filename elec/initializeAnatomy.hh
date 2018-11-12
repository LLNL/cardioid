#ifndef INITIALIZE_ANATOMY_HH
#define INITIALIZE_ANATOMY_HH

#include <mpi.h>
#include <string>

class Anatomy;

void initializeAnatomy(Anatomy& anatomy, const std::string& name, MPI_Comm comm);

#endif
