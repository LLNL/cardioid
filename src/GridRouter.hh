#ifndef GRIDROUTER_H
#define GRIDROUTER_H

#include <mpi.h>
#include <vector>
#include "Long64.hh"

class CommTable;

class GridRouter
{
 private:
   
   MPI_Comm comm_;
   int nSend_;
   std::vector<int> sendRank_;
   std::vector<int> sendOffset_;
   std::vector<int> sendIndex_;
   std::vector<std::vector<int> > instencil_;
   
  public:
  
   GridRouter(std::vector<Long64>& gid, int nx, int ny, int nz, MPI_Comm comm);
   CommTable commTable() const;
   const std::vector<int>& sendMap() const {return sendIndex_;}
};
#endif

