#ifndef GDLOADBALANCER_H
#define GDLOADBALANCER_H

#include <vector>
#include <mpi.h> 
using std::vector;

class AnatomyCell;
class ProcBox;
class MoveInfo;

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
  int gapthresh_;                // number of consecutive zero cells that defines a gap
  vector<vector<int> > penbr_;   // process numbers of all neighbors
  vector<int> thatn_;            // convenience function for neighbor pair indices
  vector<int> haveData_;         // (testing) current process owns data for processor ip
  vector<int> myRankList_;       // (testing) list of processes owned by myRank_
  vector<vector<int> > togive_;  // list of all data exchanges needed for balance
  vector<int> histnloc_;         // copy of nloc_ used for histogramming
  vector<int> histvol_;          // copy of vol_ used for histogramming
  vector<ProcBox*> peboxinfo_;  // vector of all processors' bounding box info
  vector<int> pevol_;            // bounding box volume of each processor's cells
  vector<vector<double> > pecom_;
  vector<vector<int> > pemin_;
  vector<vector<int> > pemax_;
  MPI_Comm comm_;
  int nTasks_, myRank_;
    
  public:

  GDLoadBalancer(MPI_Comm comm, int npex, int npey, int npez);
  ~GDLoadBalancer();
  void initialDistByVol(vector<AnatomyCell>& cells, int nx, int ny, int nz);
  void xstripDist(vector<int>& xcnt, vector<int>& pexmin, vector<int>& pexmax, bool verbose);
  void initialDistribution(vector<AnatomyCell>& cells, int nx, int ny, int nz);
  void gridVolMinLoop(vector<AnatomyCell>& cells);
  void assignCells(vector<AnatomyCell>& cells);
  void applyStoredMoves(vector<MoveInfo> moves,vector<AnatomyCell>& cells);
  void compactLoop(vector<AnatomyCell>& cells);
  void computeLoadInfo(vector<AnatomyCell>& cells);
  double pbcDist(double x1,double x2,int nx);
  void balanceLoop(vector<AnatomyCell>& cells);
  void balanceLoop(vector<AnatomyCell>& cells, int bblock, int bthresh, int maxiter);
  bool diffuseDest(vector<AnatomyCell>& cells);
  void redistributeCells(vector<AnatomyCell>& cells);
  void restrictMoves(vector<AnatomyCell>& cells, int ngap);
  void nonzeroVolume(vector<AnatomyCell>& cells);
  void nlocHistogram(vector<int>& nloc);
  void nlocHistogram(vector<AnatomyCell>& cells);
  void computeNlocHistogram();
  void volHistogram(vector<AnatomyCell>& cells);
  void computeVolHistogram();
  void computeVolHistogram(vector<AnatomyCell>& cells);
  void computeReducedProcGrid(int nTasks);
  void setReducedProcGrid(int rnx, int rny, int rnz);
  int ipow(int a, int b);
};

class GapSort
{
 public:
    GapSort(vector<int>& v): v_(v) {};
    vector<int> v_;
    bool operator()(const int i, const int j)
   {
      return v_[i] < v_[j];
   }
};
    
class MoveInfo
{
 public:
    MoveInfo(int spe, int dpe, int d, int v): srcpe(spe),destpe(dpe),dim(d),val(v) {};
    int srcpe;
    int destpe;
    int dim;
    int val;
};
    
#endif
