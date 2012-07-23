#include "workBound.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

double costFcn(int nTissue,int height, int dx, int dy, double a)
{
      int dz4 = (height+2)   ; 
      if (dz4 % 4  != 0) dz4 += (4-dz4%4); 
      int bbVol = (dx+2)*(dy+2)*dz4;
      double cost = nTissue+a*bbVol; 
      return cost; 
}
int  findChunks(BALANCER *balancer, unsigned short *seg, double diffCost, COLUMN *column)
{
   int NX = balancer->NX; 
   int NY = balancer->NY; 
   int nx = balancer->nx; 
   int ny = balancer->ny; 
   int nz = balancer->nz; 
   int dx = balancer->dx; 
   int dy = balancer->dy; 
   int dz = balancer->dz; 
   CHUNKS *chunkBuffer = balancer->chunkBuffer; 
  int dZ = dz+2; 
  if (dZ % 4  != 0) dZ += (4-dZ%4); 
  int maxCost = costFcn(dx*dy*dz,dZ-2,dx,dy,diffCost); 
  int offset=0; 
  int sumChunks=0; 
  for (int j=0;j<NY;j++) 
  for (int i=0;i<NX;i++) 
  {
    int sum=0; 
    int cxy = i + NX*j; 
    unsigned short *col = seg+nz*cxy; 
    int zMin ;
    int zMax ;
    for (zMin = 0   ;zMin<nz;zMin++) if (col[zMin] > 0) break ; 
    for (zMax = nz-1;zMax>=0;zMax--) if (col[zMax] > 0) break ; 
    int z0 ;
    int nChunks=0; 
    CHUNKS *chunks ; 
    if (chunkBuffer != NULL) chunks = chunkBuffer+offset;  else chunks=NULL; 
    int offsetStart = offset; 
    for (int z=zMin;z<=zMax;z++) 
    {
       int cnt = col[z]; 
       if (sum == 0) z0 = z; 
       sum+=cnt; 
       int cost = costFcn(sum,z-z0+1,dx,dy,diffCost); 
 
       if (sum > 0) 
       {
          if ( cost > maxCost  ) 
          {
              sum-=cnt; 
              cost = costFcn(sum,z-z0,dx,dy,diffCost); 
              if (chunks != NULL) 
              {
              chunks[nChunks].zL=z0; 
              chunks[nChunks].zU=z-1; 
              chunks[nChunks].nTissue=sum; 
              }
              offset++;
              nChunks++; 
              sum = cnt; 
              z0 = z; 
          }
       }
     }
     if (sum > 0) 
     {
       int dz4 = (zMax-z0+1+2)   ; 
       if (dz4 % 4  != 0) dz4 += (4-dz4%4); 
       int cost = costFcn(sum,dz4-2,dx,dy,diffCost); 
       if (chunks != NULL) 
       {
          chunks[nChunks].zL=z0; 
          chunks[nChunks].zU=zMax; 
          chunks[nChunks].nTissue=sum; 
       }
       offset++;
       nChunks++; 
     }
     if (column != NULL) 
     {
       column[cxy].nChunks = nChunks; 
       if (nChunks > 0) column[cxy].offset =  offsetStart; 
       else column[cxy].offset=-1; 
       sumChunks += nChunks; 
     }
   }
  // printf("sumChunks=%d\n",sumChunks); 
   return offset; 
}
BALANCER buildBalancer(int nx,int ny,int nz,int dx,int dy,int dz, int  nTasks, int printStats, MPI_Comm comm) 
{
   int rank,nRank; 
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &nRank);
   BALANCER balancer; 
   balancer.comm = comm; 
   balancer.printStats=printStats; 
   balancer.rank   = rank; 
   balancer.nRank  = nRank; 
   balancer.nTasks = nTasks; 
   balancer.nx = nx; 
   balancer.ny = ny; 
   balancer.nz = nz; 
   balancer.dx = dx; 
   balancer.dy = dy; 
   balancer.dz = dz; 
   balancer.NX = (nx%dx == 0) ? nx/dx : nx/dx + 1; 
   balancer.NY = (ny%dy == 0) ? ny/dy : ny/dy + 1; 
   balancer.NXYz= balancer.NX*balancer.NY*nz; 
   balancer.columns = (COLUMN *)malloc(sizeof(COLUMN)*balancer.NX*balancer.NY); 
   balancer.chunkBuffer = (CHUNKS *)malloc(sizeof(CHUNKS)*nTasks); 
   balancer.segLocal = (unsigned short *)malloc(sizeof(short unsigned)*balancer.NXYz); 
   balancer.seg      = (unsigned short *)malloc(sizeof(short unsigned)*balancer.NXYz); 
   for (int i=0;i<balancer.NXYz;i++) balancer.segLocal[i]=0; 
   return balancer; 
}
void destroyBalancer(BALANCER *balancer)
{
   free(balancer->segLocal); 
   free(balancer->seg); 
   free(balancer->chunkBuffer); 
   free(balancer->columns); 
}
void fillSeg(int nn, long long unsigned *gidArray, BALANCER *balancer) 
{
   unsigned short *seg = balancer->segLocal; 
   int rank = balancer->rank; 
   int NX = balancer->NX; 
   int NY = balancer->NY; 
   int NXYz =balancer->NXYz;
   int nx = balancer->nx; 
   int ny = balancer->ny; 
   int nz = balancer->nz; 
   int dx = balancer->dx; 
   int dy = balancer->dy; 
   int dz = balancer->dz; 
   for (int ii=0;ii<nn;ii++) 
   {
      long long unsigned  gid = gidArray[ii]; 
      int x=gid%nx; 
      int y=((gid-x)/nx)%ny; 
      int z=((gid-x-nx*y))/(nx*ny); 
      int i = x/dx;
      int j = y/dy;
      assert(0<=i && i<NX); 
      assert(0<=j && j<NY); 
      int segIndex = z + nz*(i + NX*j); 
      assert( 0 <= segIndex && segIndex < NXYz); 
      seg[segIndex] += 1; 
   }
}
void reduceSeg(BALANCER *balancer) 
{
   int NXYz=balancer->NXYz; 
   unsigned short *segSum = (unsigned short *)malloc(sizeof(unsigned short)*NXYz); 
   MPI_Allreduce(balancer->segLocal,balancer->seg,NXYz,MPI_UNSIGNED_SHORT,MPI_SUM,balancer->comm); 
}
int buildColumns(BALANCER *balancer)
{
   unsigned short *seg =balancer->seg; 
   COLUMN *columns = balancer->columns; 
   int rank=balancer->rank; 
   int NX = balancer->NX; 
   int NY = balancer->NY; 
   int nx = balancer->nx; 
   int ny = balancer->ny; 
   int nz = balancer->nz; 
   int dx = balancer->dx; 
   int dy = balancer->dy; 
   int dz = balancer->dz; 
   double minCost = -1; 
   double maxCost = -1; 
   double  relCost ; 
   int nTasks; 
   for (relCost = 0.0125;relCost<2;relCost*=2.0) 
   {
      nTasks=findChunks(balancer,seg, relCost,NULL);
      if ( nTasks < balancer->nTasks) minCost = relCost; 
      if ( nTasks > balancer->nTasks) {maxCost = relCost; break;}
   }
   if (minCost == -1) return 2;
   if (maxCost == -1) return 3; 
   int loop=0; 
   while(nTasks != balancer->nTasks   )
   {
     relCost = 0.5*(minCost+maxCost); 
     nTasks=findChunks(balancer,seg, relCost, NULL);
     if (nTasks < balancer->nTasks) minCost = relCost; 
     if (nTasks > balancer->nTasks) maxCost = relCost; 
     if (loop++ > 30) break; 
    }
    balancer->relCost = relCost;
   if (nTasks == balancer->nTasks)  
   {
      nTasks=findChunks(balancer,seg, relCost, columns);
      return  0;
   }
   return 1; 
}
void  printDomainInfo(BALANCER *balancer) 
{
   if (balancer->rank != 0) return ; 
   FILE *file=fopen("domainInfo.data","w"); 
   COLUMN *column=balancer->columns; 
   for (int cxy =0;cxy<balancer->NX*balancer->NY;cxy++) 
   {
     int offset = column[cxy].offset; 
     for (int ic=0;ic<column[cxy].nChunks;ic++) 
     {
     
       CHUNKS* chunks = balancer->chunkBuffer+(offset);
        
       int height = chunks[ic].zU-chunks[ic].zL+1;
       int dz4 = (height+2)   ; 
       if (dz4 % 4  != 0) dz4 += (4-dz4%4); 
       int bbVol = (balancer->dx+2)*(balancer->dy+2)*dz4;
       int nTissue = chunks[ic].nTissue; 
       double cost = nTissue+balancer->relCost*bbVol; 
       int id= offset+ic; 
       fprintf(file,"%8d %4d %6d %4d %4d %f\n",id,nTissue,bbVol,height,dz4,cost); 
     }
   }
   fclose(file); 
}
BALANCE_DOMAIN  domainInfo(long long unsigned  gid, BALANCER *balancer) 
{
   COLUMN *column=balancer->columns; 
   int x=gid%balancer->nx; 
   int y=((gid-x)/balancer->nx)%balancer->ny; 
   int z=((gid-x-balancer->nx*y))/(balancer->nx*balancer->ny); 
   int i = x/balancer->dx;
   int j = y/balancer->dy;
   int cxy = (i + balancer->NX*j); 
   BALANCE_DOMAIN domain; 
   for (int ic=0;ic<column[cxy].nChunks;ic++) 
   {
       CHUNKS* chunks = balancer->chunkBuffer+(column[cxy].offset);
      if (chunks[ic].zL <= z &&  z<= chunks[ic].zU) 
      {
          int height = chunks[ic].zU-chunks[ic].zL+1;
          int dz4 = (height+2)   ; 
          if (dz4 % 4  != 0) dz4 += (4-dz4%4); 
          int bbVol = (balancer->dx+2)*(balancer->dy+2)*dz4;
          domain.bbVol =bbVol; 
          domain.bbHeight = dz4; 
          domain.id= column[cxy].offset+ic; 
          return domain; 
      }
   }
  domain.bbHeight =0; 
  domain.bbVol =0; 
  domain.id= -1  ;
  return domain; 
}
