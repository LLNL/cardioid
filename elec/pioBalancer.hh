#ifndef PIO_BALANCER_HH
#define PIO_BALANCER_HH

#include <string>
#include <mpi.h>

class Simulate;

int pioBalancer(const std::string& domainFile,
                const std::string& pxyzFile,
                Simulate& sim,
                MPI_Comm comm);
#endif 
