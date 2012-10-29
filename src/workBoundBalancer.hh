#ifndef WORKBOUNDBALANCER_HH
#define WORKBOUNDBALANCER_HH
#include "workBound.h"
#include "LoadLevel.hh"
LoadLevel workBoundBalancer(vector<AnatomyCell> &cells, int dx, int dy, int dz, int nx, int ny, int nz, int target, 
                      int nCores, int nRCoresBB, double alpha, double beta, int printStats,MPI_Comm comm);
#endif
