#include "workBound.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>
#define MIN(A,B) ((A) < (B) ? (A) : (B))

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
    if (column != NULL) chunks = chunkBuffer+offset;  else chunks=NULL; 
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
              chunks[nChunks].domainID=offset; 
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
          chunks[nChunks].domainID=offset; 
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
   return offset; 
}
BALANCER buildBalancer(long long unsigned nCellGlobal, int nx,int ny,int nz,int dx,int dy,int dz, 
       int  nTasks, int nCores, int nRCoresBB, double alpha, double beta, int printStats, MPI_Comm comm) 
{
   int rank,nRank; 
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &nRank);
   int dz4 = (dz+2)   ;
   if (dz4 % 4  != 0) dz4 += (4-dz4%4);

   BALANCER balancer; 
   balancer.comm = comm; 
   balancer.nCellGlobal = nCellGlobal; 
   balancer.printStats=printStats; 
   balancer.nC = nCores; 
   balancer.alphaWork=alpha; 
   balancer.beta=beta; 
   //balancer.timeRef = (dx*dy*dz + alpha*(dx)*(dy)*(dz4))/nC ;
   balancer.timeRef = (1.0*dx*dy*dz)/(nRCoresBB) ;
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
   unsigned short *segLocal = balancer->segLocal; 
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
      segLocal[segIndex] += 1; 
   }
}
void reduceSeg(BALANCER *balancer) 
{
   int NXYz=balancer->NXYz; 
   MPI_Allreduce(balancer->segLocal,balancer->seg,NXYz,MPI_UNSIGNED_SHORT,MPI_SUM,balancer->comm); 
}
void  setDomainInfo(BALANCER *balancer) 
{
   COLUMN *column=balancer->columns; 
   for (int cxy =0;cxy<balancer->NX*balancer->NY;cxy++) 
   {
     int offset = column[cxy].offset; 
     CHUNKS* chunks = balancer->chunkBuffer+(offset);
     for (int ic=0;ic<column[cxy].nChunks;ic++) 
     {
       int  id= offset+ic; 
       if (id == balancer->rank)
       {
           balancer->chunk = chunks[ic];
           break; 
       }
     }
   }
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
   double  alpha ; 
   int nTasks; 
   for (alpha = 0.0125;alpha<2;alpha*=2.0) 
   {
      nTasks=findChunks(balancer,seg, alpha,NULL);
      if (nTasks == balancer->nTasks)  {minCost=maxCost=alpha; break;}
      if ( nTasks < balancer->nTasks) minCost = alpha; 
      if ( nTasks > balancer->nTasks) {maxCost = alpha; break;}
   }
   //if (minCost == -1) return 2;
   //if (maxCost == -1) return 3; 
   int loop=0; 
   while(nTasks != balancer->nTasks   )
   {
     alpha = 0.5*(minCost+maxCost); 
     nTasks=findChunks(balancer,seg, alpha, NULL);
     if (nTasks < balancer->nTasks) minCost = alpha; 
     if (nTasks > balancer->nTasks) maxCost = alpha; 
     if (loop++ > 30) break; 
    }
    balancer->alphaBalance = alpha;
   if (nTasks == balancer->nTasks)  
   {
      nTasks=findChunks(balancer,seg, alpha, columns);
      return  0;
   }
   return 1; 
}
void  printDomainInfo(BALANCER *balancer) 
{
   if (balancer->rank != 0) return ; 
   char filename[256]; 
   double nAve = balancer->nCellGlobal/balancer->nTasks; 
   double nMax = balancer->dx*balancer->dy*balancer->dz; 
   sprintf(filename,"domains#%6.6d",balancer->rank);
   FILE *domainFile = fopen(filename,"w");
   fprintf(domainFile,"domain HEADER {\n datatype = VARRECORDASCII;\n");
   fprintf(domainFile,"nfiles=1; nrecords=%d; nfields=15;\n",balancer->nTasks);
   fprintf(domainFile,"nfields=15;\n");
   fprintf(domainFile,"field_names=domain cx cy zL ZU nT nBB, dz dz4 cost nR nD timeR timeD timeRef;\n" );
   fprintf(domainFile,"field_types=u u u u u u u u u f u u f f f;\n" );
   fprintf(domainFile,"nCells=%llu; nx=%d; ny=%d; nz=%d;\n",balancer->nCellGlobal,balancer->nx,balancer->ny,balancer->nz);
   fprintf(domainFile,"dx=%d; dy=%d; dz=%d; ",balancer->dx,balancer->dy,balancer->dz);
   fprintf(domainFile,"NX=%d; NY=%d;\n",balancer->NX,balancer->NY);
   fprintf(domainFile,"alphaWork=%f; alphaBalance=%f; beta=%f\n",balancer->alphaWork,balancer->alphaBalance,balancer->beta);
   fprintf(domainFile,"nCellAve=%f; nCellBox=%f; ratio=%f;\n", nAve,nMax,nMax/nAve);
   fprintf(domainFile,"nCores=%d;\n}\n\n",balancer->nC);
 
   int cnt=0;
   double timeRef = balancer->timeRef; 
   for (int cxy =0;cxy<balancer->NX*balancer->NY;cxy++)
   {
      int x = cxy % balancer->NX;
      int y = cxy / balancer->NX;
      CHUNKS *chunks = (balancer->chunkBuffer+(balancer->columns[cxy].offset));
      for (int j=0;j<balancer->columns[cxy].nChunks;j++)
      {
        CHUNKS chunk = chunks[j];
        int zL = chunk.zL;
        int zU = chunk.zU;

        PARTITION_INFO part = corePartition(chunk, balancer);
        int nT = part.nT;
        int nB = part.nB;
        int height = part.height;
        int dz4   = part.dz4;
        int nR = part.nR;
        int nD = part.nD;
        double cost = part.cost;
        double timeR = part.timeR;
        double timeD = part.timeD;

        fprintf(domainFile,"%6d %6d %6d %6d %6d %6d %6d %6d %6d %9.3f %6d %6d %9.3f %9.3f %9.3f ",
          chunk.domainID,x,y,zL,zU,nT,nB,height,dz4,cost,nR,nD,timeR,timeD,timeRef);
        if ( timeR > 1.001*timeRef || timeD > 1.001*timeRef) fprintf(domainFile," *"); 
        fprintf(domainFile,"\n"); 

        cnt ++;

      }
   }
   fclose(domainFile); 
}
int  domainID(int gid, BALANCER *balancer)
{
   COLUMN *column = balancer->columns;
   int x=gid%balancer->nx;
   int y=((gid-x)/balancer->nx)%balancer->ny;
   int z=((gid-x-balancer->nx*y))/(balancer->nx*balancer->ny);
   int i = x/balancer->dx;
   int j = y/balancer->dy;
   int cxy = (i + balancer->NX*j);
   for (int ic=0;ic<column[cxy].nChunks;ic++)
   {
       CHUNKS* chunks = balancer->chunkBuffer+(column[cxy].offset);
      if (chunks[ic].zL <= z &&  z<= chunks[ic].zU) return chunks[ic].domainID;
   }
   assert(0);
   return -1;
}

int balanceCores(BALANCER *balancer)
{
   CHUNKS chunk = balancer->chunk;
   int nCore = balancer->nC; 
   int dx = balancer->dx;
   int dy = balancer->dy;
   int dz = balancer->dz;
   int nTissue = chunk.nTissue; 
   int height =  chunk.zU-chunk.zL+1;
   int dz4 = (height+2)   ;
   if (dz4 % 4  != 0) dz4 += (4-dz4%4);
   int bbVol = (dx+2)*(dy+2)*dz4;
   int nDiffusion;
   int hh= dz4 /4;
   if  (hh <=4 ) nDiffusion = 2;
   if  (hh > 4 ) nDiffusion = hh/2;
   int nReaction = nCore-nDiffusion;

   int volR = dx*dy*dz;
   int volD = (dx+2)*(dy+2)*(dz*2); 
   double timeR = (1.0*nTissue)/volR;
   double timeD = (1.0/7.0*bbVol)/volR * (nCore-nDiffusion)/nDiffusion; 
   while (timeD < 0.95 && timeR > 1.05 && nDiffusion > 1)
   {
      nReaction++;
      nDiffusion--;
      timeR = (1.0*nTissue)/volR;
      timeD = (1.0/7.0*bbVol)/volR * (nCore-nDiffusion)/nDiffusion; 
   }
   while (timeR < 0.95 && timeD > 1.05 && nReaction > 1)
   {
      nReaction--;
      nDiffusion++;
      timeR = (1.0*nTissue)/volR;
      timeD = (1.0/7.0*bbVol)/volR * (nCore-nDiffusion)/nDiffusion; 
   }
    return nDiffusion; 
}
PARTITION_INFO corePartition(CHUNKS chunk, BALANCER *balancer)
{
   PARTITION_INFO part; 
   int dx = balancer->dx;
   int dy = balancer->dy;
   int dz = balancer->dz;
   int nC = balancer->nC;
   double alpha = balancer->alphaWork;
   double beta = balancer->beta;
   double timeRef = balancer->timeRef;
   part.nT = chunk.nTissue;
   part.height =  chunk.zU-chunk.zL+1;
   part.dz4 = (part.height+2)   ;
   if (part.dz4 % 4  != 0) part.dz4 += (4-part.dz4%4);
   part.nB = (dx)*(dy)*part.dz4;
   double timeMin =  20*part.nT;
   int nRMin=1;
   for (int nR =1;nR<nC;nR++)
   {
       int nD = nC-nR;
       double timeD_simd    = (alpha*part.nB)/nD;
       double timeD_threads = (alpha*beta*part.nT)/nD;
       double timeD = MIN(timeD_simd,timeD_threads); 
       double timeR = (1.0*part.nT)/nR;
       double time = (timeD > timeR) ?  timeD : timeR;
       if ( time < timeMin ) {timeMin=time;nRMin=nR;}
   }
   int nR = nRMin; 
   int nD = nC-nR; 
   part.nR = nR;
   part.nD = nD;
   double timeD_simd    = (alpha*part.nB)/nD;
   double timeD_threads = (alpha*beta*part.nT)/nD;
   double timeD = MIN(timeD_simd,timeD_threads); 
   
   if (timeD_simd > timeD_threads) part.stencilMethod = strdup("strip"); 
   else                            part.stencilMethod = strdup("threads"); 
   part.timeD = timeD;
   part.timeR = (1.0  *part.nT)/part.nR;
   part.cost  = part.nT + alpha*part.nB;
   return part; 
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
