#ifndef ASSIGN_CELLS_TO_TASKS_HH
#define ASSIGN_CELLS_TO_TASKS_HH

#include <string>
#include <mpi.h>
#include "loadlevel.hh" 

class Simulate;

LoadLevel assignCellsToTasks(Simulate& sim, const std::string& name, MPI_Comm comm);

#endif
