#include <mpi.h>
#include <stdio.h>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <vector>

#include "mpiUtils.h"
#include "Simulate.hh"
#include "AnatomyCell.hh"
#include "workBound.h"
using namespace std;
int workBoundAssignCells(vector<AnatomyCell> &cells, int dx, int dy, int dz, int nx, int ny, int nz, int target)
{
   char filename[16]; 

   BALANCER balancer =buildGeom(nx,ny,nz,dx,dy,dz,target);
   unsigned nLocal = cells.size(); 
   int *gid= (int *)malloc(sizeof(int)*nLocal); 
   for (int i=0;i<nLocal;i++) gid[i]=cells[i].gid_; 
   unsigned short *seg=buildSeg(nLocal,gid,&balancer);
   //sprintf(filename,"Stuff#%6.6d",balancer.rank);
   //FILE *stuff = fopen(filename,"w");
   COLUMN *columns = buildColumns(seg,&balancer);
   destroySeg(seg);
   for (int i =0;i<nLocal;i++) 
   {
     int id = domainID(gid[i], columns, &balancer) ;
     cells[i].dest_= id ;
    // fprintf(stuff,"%d %d\n",gid[i],id);
   }
   free(gid); 
   unsigned nGlobal; 
   MPI_Reduce(&nLocal,&nGlobal,1,MPI_UNSIGNED,MPI_SUM,0,MPI_COMM_WORLD); 
   if (balancer.rank ==0) printf("nGlobal=%d\n",nGlobal); 
   MPI_Barrier(MPI_COMM_WORLD); 


   AnatomyCellDestSort sortByDest;
   sort(cells.begin(), cells.end(), sortByDest);
   vector<unsigned> dest(nLocal);
   for (int i =0;i<nLocal;i++) dest[i] = cells[i].dest_ ;//domainID( (int)cells[i].gid_, columns, &balancer) ;
   destroyColumns(columns);

   unsigned capacity = max((int)nLocal,dx*dy*dz);
   cells.resize(capacity);
   assignArray((unsigned char*) &(cells[0]),
               &nLocal,
               capacity,
               sizeof(AnatomyCell),
               &(dest[0]),
               0,
               MPI_COMM_WORLD);
   cells.resize(nLocal);
   MPI_Reduce(&nLocal,&nGlobal,1,MPI_UNSIGNED,MPI_SUM,0,MPI_COMM_WORLD); 
   if (balancer.rank ==0) printf("nGlobal=%d\n",nGlobal); 
   sort(cells.begin(), cells.end(), sortByDest);
   int nDiffusionTasks = -1; 
   return nDiffusionTasks; 
}

