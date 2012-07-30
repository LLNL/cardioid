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
#include "pio.h"
using namespace std;
int workBoundBalancer(vector<AnatomyCell> &cells, int dx, int dy, int dz, int nx, int ny, int nz, int target,int nC, double alpha, int printStats,MPI_Comm comm)
{
   unsigned nLocal = cells.size(); 
   unsigned nGlobal; 
   MPI_Allreduce(&nLocal,&nGlobal,1,MPI_UNSIGNED,MPI_SUM,comm); 
   BALANCER balancer = buildBalancer(nGlobal,nx,ny,nz,dx,dy,dz,target,nC,alpha,printStats,comm);

   int rank = balancer.rank; 
   //char filename[16]; 
   //sprintf(filename,"Stuff#%6.6d",rank);
   //FILE *stuff = fopen(filename,"w");

   unsigned int nn = 1024; 
   vector<long long unsigned> gid(nn); 
   for (unsigned int ii=0;ii<nLocal;ii+=nn) 
   {
       unsigned int n = ((nLocal-ii) >= nn) ? nn : (nLocal-ii) ; 
       for (unsigned int j=0;j<n;j++) gid[j]=(cells[ii+j].gid_); 
       fillSeg(n,&gid[0],&balancer);
   }
   reduceSeg(&balancer);

   int rc = buildColumns(&balancer);
   if (rank==0 && rc == 1) printf("buildColumns failed to converge to target number of tasks"); 
   if (rank==0 && rc == 2) printf("buildColumns failed to minA"); 
   if (rank==0 && rc == 3) printf("buildColumns failed to maxA"); 
   if (rc !=0) exit(0); 

   for (int i =0;i<nLocal;i++) 
   {
      BALANCE_DOMAIN domain = domainInfo(cells[i].gid_, &balancer) ;
      cells[i].dest_ = domain.id ;
   }
   if (balancer.printStats)
   {
      PFILE* file = Popen("domainMap.data", "w", balancer.comm);
      if (balancer.rank ==0) Pprintf(file,"domainMap HEADER {labels=gid task}\n\n"); 
      for (int i =0;i<nLocal;i++) Pprintf(file,"%llu %u\n",cells[i].gid_,cells[i].dest_); 
      Pclose(file); 
      printDomainInfo(&balancer); 
   }

   if (balancer.nRank != balancer.nTasks) 
   {
       if ( balancer.rank == 0) printf("target number of tasks differ from actual number. Exiting Program\n"); 
         exit(0);   
   }
   int nDiffusionTasks = -1; 
   setDomainInfo(&balancer); 
   PARTITION_INFO part = corePartition(balancer.chunk,&balancer);
   nDiffusionTasks = part.nD; 
   destroyBalancer(&balancer); 

   MPI_Allreduce(&nLocal,&nGlobal,1,MPI_UNSIGNED,MPI_SUM,balancer.comm); 
   if (rank==0) {printf("Before particle assignment nGlobal=%d\n",nGlobal); fflush(stdout);}
   MPI_Barrier(balancer.comm); 

   AnatomyCellDestSort sortByDest;
   sort(cells.begin(), cells.end(), sortByDest);
   vector<unsigned> dest(nLocal);
   for (int i =0;i<nLocal;i++) dest[i] = cells[i].dest_ ;

   unsigned capacity = max((int)nLocal,dx*dy*dz);
   cells.resize(capacity);
   assignArray((unsigned char*) &(cells[0]), &nLocal, capacity,
               sizeof(AnatomyCell), &(dest[0]), 0, balancer.comm);
   cells.resize(nLocal);
   sort(cells.begin(), cells.end(), sortByDest);

   MPI_Allreduce(&nLocal,&nGlobal,1,MPI_UNSIGNED,MPI_SUM,balancer.comm); 
   if (rank==0) {printf("After particle assignment nGlobal=%d\n",nGlobal); fflush(stdout);}

   return nDiffusionTasks; 
}

