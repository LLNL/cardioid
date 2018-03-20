#ifndef TAGSERVER_H
#define TAGSERVER_H

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Returns a unique tag that can be used for MPI communication.*/
int getUniqueTag(MPI_Comm comm);

#ifdef __cplusplus
}
#endif

#endif
