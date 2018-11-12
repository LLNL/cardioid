#ifndef GRIDROUTER_H
#define GRIDROUTER_H

#include <mpi.h>
#include <vector>
#include "Long64.hh"

class CommTable;

class GridRouter
{
 private:

  int selfTest();
   
   MPI_Comm comm_;
   std::vector<int> sendRank_;
   std::vector<int> sendOffset_;
   std::vector<int> sendMap_;
   
  public:
  
   GridRouter(std::vector<Long64>& gid, int nx, int ny, int nz, MPI_Comm comm);
   CommTable commTable() const;
   const std::vector<int>& sendMap() const {return sendMap_;}
};
#endif

