#ifndef  WORKBOUNDBALANCER_HH
#define WORKBOUNDBALANCER_HH
#include "workBound.h"
int workBoundBalancer(vector<AnatomyCell> &cells, int dx, int dy, int dz, int nx, int ny, int nz, int target, int nC, double alpha, int printStats,MPI_Comm comm);
#endif
