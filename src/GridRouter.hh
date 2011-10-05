#ifndef GRIDROUTER_H
#define GRIDROUTER_H

#include <mpi.h>
#include <vector>
using namespace std;

class CommTable;

class GridRouter
{
  private:

  MPI_Comm comm_;
  int nSend_;
  vector<int> sendRank_;
  vector<int> sendOffset_;
  vector<int> sendIndex_;
  vector<vector<int> > instencil_;
  
  public:
  
  GridRouter(vector<int>& locgid, int nx, int ny, int nz, vector<double>& center, double radius, MPI_Comm comm);
  ~GridRouter(void);
  CommTable* fillTable(void);
  
};
#endif

