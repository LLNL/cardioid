#ifndef GDLOADBALANCER_H
#define GDLOADBALANCER_H

#include <vector>
#include <mpi.h> 


class GDLoadBalancer {
  private:

  int ntissue_;                  // number of non-zero grid points
  int npex_, npey_, npez_;       // dimensions of process grid
  int nx_, ny_, nz_;             // dimensions of data grid
  int npegrid_;                  // size of process grid
  int nnbr_;                     // number of neighboring processors to share data with
  int nloctot_;
  double nlocavg_;
  std::vector<int> nloc_;                  // number of non-zero grid points on each process
  std::vector<std::vector<int> > penbr_;   // process numbers of all neighbors
  std::vector<int> thatn_;                 // convenience function for neighbor pair indices
  std::vector<int> gpe_;                   // process that owns each grid point (may want to ditch this)
  std::vector<std::vector<int> > togive_;  // list of all data exchanges needed for balance
  std::vector<std::vector<int> > loctype_; // output: type list for each pe's local data
  std::vector<std::vector<int> > locgid_;  // output: gid list for each pe's local data
  
  public:

  GDLoadBalancer(int npex, int npey, int npez);
  ~GDLoadBalancer();
  void initialDistribution(std::vector<int>& types, int nx, int ny, int nz);
  void balanceLoop(void);
  void balanceLoop(int bblock, int bthresh, int maxiter);
  const std::vector<std::vector<int> >& loctype(void) const { return loctype_; }
  const std::vector<std::vector<int> >& locgid(void) const { return locgid_; }
  void loadHistogram();
  
};
#endif
