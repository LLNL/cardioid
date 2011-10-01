#ifndef GDLOADBALANCER_H
#define GDLOADBALANCER_H

#include <vector>
#ifdef USE_MPI
#include <mpi.h> 
#else
typedef int MPI_Comm;
#endif
using namespace std;

class GDLoadBalancer {
  private:

  int ntissue_;                  // number of non-zero grid points
  int npex_, npey_, npez_;       // dimensions of process grid
  int npegrid_;                  // size of process grid
  int nnbr_;                     // number of neighboring processors to share data with
  int nloctot_;
  double nlocavg_;
  vector<int> nloc_;             // number of non-zero grid points on each process
  vector<vector<int> > penbr_;   // process numbers of all neighbors
  vector<int> thatn_;            // convenience function for neighbor pair indices
  vector<int> gpe_;              // process that owns each grid point (may want to ditch this)
  vector<vector<int> > loctype_; // type list for each pe's local data
  vector<vector<int> > locgid_;  // gid list for each pe's local data
  vector<vector<int> > tosend_;  // output: list of all data exchanges needed for balance
  
  public:

  GDLoadBalancer(int npex, int npey, int npez);
  ~GDLoadBalancer();
  void initialDistribution(vector<int>& types, int nx, int ny, int nz);
  void balanceLoop(void);
  void balanceLoop(int bblock, int bthresh, int maxiter);
  const std::vector<std::vector<int> >& tosend(void) const { return tosend_; }
  void loadHistogram();
  
};
#endif
