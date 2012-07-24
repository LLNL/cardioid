#ifndef BLOCKLOADBALANCER_H
#define BLOCKLOADBALANCER_H

#include <vector>
#include <mpi.h> 
using std::vector;

class AnatomyCell;

class BlockLoadBalancer {

  private:

  int nx_, ny_, nz_;             // dimensions of data grid
  int bx_, by_, bz_;             // optimal block size
  MPI_Comm comm_;
  int nTasks_, myRank_;

  public:

    BlockLoadBalancer(MPI_Comm comm, int nx, int ny, int nz, int bx, int by, int bz);
    ~BlockLoadBalancer();
    int block(vector<AnatomyCell>& cells, double diffCost, int nCols, vector<int>& myCols);
    double costFunction(int nTissue, int height, double a);
};

#endif
