#ifndef TAGSERVER_H
#define TAGSERVER_H

#include <mpi.h>

/** Returns a unique tag that can be used for MPI communication.*/
int getUniqueTag(MPI_Comm comm);

#endif
