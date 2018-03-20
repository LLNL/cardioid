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
// The next four functions return meaningful values only when
// hi_hasTorus is true.
int  hi_nTorusDim(void);
// caller must ensure arrays are large enough to store nTorusDim
// items
void hi_torusCoords(int* coord); 
void hi_torusSize(int* size);
int hi_coreId(void);

#ifdef __cplusplus
}
#endif

#endif // #ifndef HARDWAREINFO_H


/* Local Variables: */
/* tab-width: 3 */
/* End: */
