#include "readSnapshotCellList.hh"
#include "Simulate.hh"
#include "mpi.h"
#include <iostream>
#include <fstream>
using namespace std;

bool readSnapshotCellList(string filename, Simulate& sim)
{
   int nTasks, myRank;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   
   // try to open file
   int openfail;
   ifstream input;
   if (myRank == 0)
   {
      input.open(filename.c_str(),ifstream::in);
      openfail = 0;
      if (!input.is_open())
         openfail = 1;
   }
   MPI_Bcast(&openfail, 1, MPI_INT, 0, MPI_COMM_WORLD);
   if (openfail == 1)
   {
      if (myRank == 0)
         cout << "Could not open cell list file " << filename << endl;
      return false;
   }

   // read gids from file
   vector<Long64> cellVec_;
   int nSubset;
   if (myRank == 0)
   {
      while (!input.eof()) {
         Long64 igid;
         input >> igid;
         cellVec_.push_back(igid);
      }
      nSubset = cellVec_.size();
   }   
   MPI_Bcast(&nSubset, 1, MPI_INT, 0, MPI_COMM_WORLD);
   cellVec_.resize(nSubset);
   MPI_Bcast(&cellVec_[0], nSubset, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

   // now every task knows the full list of cells, want to populate sim.snapshotCellList
   // with the cells local to each task.  Copy local cell gids into an STL::Set for faster
   // searching
   set<Long64> cellSet_;
   Long64 nLocal = sim.anatomy_.nLocal();
   const vector<AnatomyCell>& cells = sim.anatomy_.cellArray();
   for (unsigned ii=0; ii<nLocal; ++ii)
   {
      Long64 gid = cells[ii].gid_;
      cellSet_.insert(gid);
   }

   for (unsigned jj=0; jj<nSubset; ++jj)
   {
      Long64 gid = cellVec_[jj];
      set<Long64>::iterator it = cellSet_.find(gid);
      if (it != cellSet_.end())
         sim.snapshotCellList_.insert(gid);
   }
   
   return true;
}
