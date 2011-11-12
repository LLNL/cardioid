/* $Id$ */ 
#ifndef HARDWAREINFO_H
#define HARDWAREINFO_H

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

int  hi_nIoTasks(MPI_Comm comm);
const int* hi_ioTaskList(MPI_Comm comm);

int hi_hasTorus(void);
void hi_torusCoords(int* x, int* y, int* z, int* t);
void hi_torusSize(int* x, int* y, int* z);
#ifdef __cplusplus
}
#endif

#endif // #ifndef HARDWAREINFO_H


/* Local Variables: */
/* tab-width: 3 */
/* End: */
