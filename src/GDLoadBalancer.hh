#ifndef GDLOADBALANCER_H
#define GDLOADBALANCER_H

#include <vector>
#include <mpi.h> 
using std::vector;

class AnatomyCell;

class GDLoadBalancer {

  private:

  int ntissue_;                  // number of non-zero grid points
  int npex_, npey_, npez_;       // dimensions of process grid
  int tnx_, tny_, tnz_;          // reduced process grid for test runs
  int nx_, ny_, nz_;             // dimensions of data grid
  int npegrid_;                  // size of process grid
  int nnbr_;                     // number of neighboring processors to share data with
  int nloctot_;
  double nlocavg_;
  vector<vector<int> > penbr_;   // process numbers of all neighbors
  vector<int> thatn_;            // convenience function for neighbor pair indices
  vector<int> haveData_;         // (testing) current process owns data for processor ip
  vector<int> myRankList_;       // (testing) list of processes owned by myRank_
  vector<vector<int> > togive_;  // list of all data exchanges needed for balance
  vector<int> histloc_;          // copy of nloc_ used for histogramming
  vector<int> pevol_;            // bounding box volume of each processor's cells
  vector<vector<double> > pecom_;
  vector<vector<int> > pemin_;
  vector<vector<int> > pemax_;
  MPI_Comm comm_;
  int nTasks_, myRank_;
    
  public:

  GDLoadBalancer(MPI_Comm comm, int npex, int npey, int npez);
  ~GDLoadBalancer();
  void initialDistribution(vector<AnatomyCell>& cells, int nx, int ny, int nz);
  void compactLoop(vector<AnatomyCell>& cells);
  void computeProcBoxes(vector<AnatomyCell>& cells);
  double pbcDist(double x1,double x2,int nx);
  void balanceLoop(vector<AnatomyCell>& cells);
  void balanceLoop(vector<AnatomyCell>& cells, int bblock, int bthresh, int maxiter);
  bool diffuseDest(vector<AnatomyCell>& cells);
  void redistributeCells(vector<AnatomyCell>& cells);
  void restrictMoves(vector<AnatomyCell>& cells, int ngap);
  void loadHistogram(vector<int>& nloc);
  void loadHistogram(vector<AnatomyCell>& cells);
  void computeHistogram();
  void computeReducedProcGrid(int nTasks);
  int ipow(int a, int b);
    
};
#endif
