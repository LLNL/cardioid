// $Id$

#ifndef MPI_UTILS_H
#define MPI_UTILS_H

#ifdef WITH_MPI
#include <mpi.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifdef WITH_MPI
void timestampBarrier(const char* msg, MPI_Comm comm);

void distributeArray(unsigned char* data,
		     const unsigned nHave,
		     const unsigned nWant,
		     const unsigned width,
		     MPI_Comm comm);

void assignArray(unsigned char* data,
		 unsigned* nLocal,
		 const unsigned capacity,
		 const unsigned width,
		 const unsigned* dest,
		 int verbose,
		 MPI_Comm comm);


#endif // #ifdef WITH_MPI

int getRank(int flag); 
int getSize(int flag); 
void WAIT(int flag);
void abortAll_core(int rc,char *file,int line);

#define abortAll(x) abortAll_core((x),__FILE__,__LINE__)

#ifdef __cplusplus
}
#endif
#endif // #ifndef MPI_UTILS_H



/* Local Variables: */
/* tab-width: 3 */
/* End: */
