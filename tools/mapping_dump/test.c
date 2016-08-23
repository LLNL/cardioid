#include <stdio.h>
#include <inttypes.h>
#include "mpi.h"

#include "spi_impl.h"

void spi_dump_mapping(spi_hdl_t* spi_hdl);
uint32_t mapping_table_new(spi_hdl_t* spi_hdl);

#define timebase(x) asm volatile ("mftb %0" : "=r"(x) : )

int main(int argc, char** argv)
{
   uint64_t t1,t2,t3,t4;
   int npes, mype;
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &npes);
   MPI_Comm_rank(MPI_COMM_WORLD, &mype);

   spi_hdl_t spi_hdl;
   mapping_table_new(&spi_hdl);
   spi_dump_mapping(&spi_hdl);

  MPI_Finalize();
  return 0;
}
