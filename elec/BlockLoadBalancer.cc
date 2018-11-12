#include "BlockLoadBalancer.hh"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <cassert>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <map>
#include <mpi.h>
#include "mpiUtils.h"
#include "GridPoint.hh"
#include "AnatomyCell.hh"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
BlockLoadBalancer::BlockLoadBalancer(MPI_Comm comm, int nx, int ny, int nz,
                                     int bx, int by, int bz): comm_(comm),
                                                              bx_(bx), by_(by), bz_(bz),
                                                              nx_(nx),ny_(ny),nz_(nz)
{
   MPI_Comm_size(comm_, &nTasks_);
   MPI_Comm_rank(comm_, &myRank_);
}
////////////////////////////////////////////////////////////////////////////////
BlockLoadBalancer::~BlockLoadBalancer()
{
}
////////////////////////////////////////////////////////////////////////////////
int BlockLoadBalancer::block(vector<AnatomyCell>& cells, double diffCost, int nCols, vector<int>& myCols)
{
   // cells should contain only the cells with non-zero type indices, distributed
   // in columns of grid points as follows:
   //   col_i = x/bx, col_j = y/by
   //   data owner pe = col_i + col_j*nblkx
   
   int nLocal = cells.size();
   int ncx = (nx_%bx_ == 0 ? nx_/bx_ : nx_/bx_ + 1);  // number of columns in x
   int ncy = (ny_%by_ == 0 ? ny_/by_ : ny_/by_ + 1);  // number of columns in y
   const double maxCost = costFunction(bx_*by_*bz_,bz_,diffCost);
   const int zlen = 4;
   
   int npes = 0;   
   if (nLocal > 0)
   {
      int pecnt = 0;
      for (int ic=0; ic<myCols.size(); ic++)
      {
         int pecol = myCols[ic];

         // sort data column by z
         int nfound = 0;
         for (unsigned ii=0; ii<cells.size(); ++ii)
         {
            GridPoint gpt(cells[ii].gid_,nx_,ny_,nz_);
            int iicol = (gpt.x/bx_) + ncx*(gpt.y/by_);
            if (iicol == pecol)
            {
               cells[ii].sortind_ = gpt.z;
               nfound++;
            }
            else
               cells[ii].sortind_ = -1;
            
         }
         if (nfound > 0)
         {
            sort(cells.begin(),cells.end(),AnatomyCell::indLessThan);

            int addcnt = 0;
            GridPoint first(cells[nLocal-nfound].gid_,nx_,ny_,nz_);
            int zmin = first.z;
            int zmax = zmin + zlen;
            for (unsigned ii=nLocal-nfound; ii<nLocal; ++ii)
            {
               addcnt++;
               GridPoint gpt(cells[ii].gid_,nx_,ny_,nz_);
               if (gpt.z < zmax)
                  cells[ii].dest_ = npes;
               else
               {
                  while (zmax <= gpt.z) zmax++;
                  double cost = costFunction(addcnt,gpt.z-zmin+1,diffCost);
                  if (cost > maxCost)
                  {
                     zmin = gpt.z;
                     zmax = zmin + zlen;
                     npes++;
                     addcnt = 1;
                  }
                  cells[ii].dest_ = npes;
               }
            }
            npes++; // subsequent columns should not share a pe with this one
         }       
      }
   }
   return npes;
}    
////////////////////////////////////////////////////////////////////////////////
double BlockLoadBalancer::costFunction(int nTissue, int height, double a)
{
   int dz4 = (height+2)   ;
   if (dz4 % 4  != 0) dz4 += (4-dz4%4);
   int bbVol = (bx_+2)*(by_+2)*dz4;
   double cost = nTissue+a*bbVol;
   return cost;
}
