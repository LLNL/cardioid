#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <unistd.h>

#include "clocksync.h"

int main(int argc,char *argv[]) {
  const int nmsg = 2000;
  int np,pid;
  int iter;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&np);
  MPI_Comm_rank(MPI_COMM_WORLD,&pid);

  cs_inittime();
  (void) cs_gettime();

  for(iter = 0; iter < 5; iter++) {
    if(pid == 0) printf("%% @@ SYNC ITERATION %d\n",iter);

    if(iter > 0)
      sleep(2);

    (void) cs_clocksync(np,pid,nmsg,
			NULL,NULL);
  }

  MPI_Finalize();

  return 0;
}
