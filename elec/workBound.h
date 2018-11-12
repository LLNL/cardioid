#ifndef WORKBOUND_H
#define WORKBOUND_H
#include <mpi.h>

typedef  struct partition_info_st { int nD, nR, nT, nB, height, dz4, cost, timeR, timeD; char *stencilMethod;} PARTITION_INFO;
typedef  struct chunks_st { int  zL, zU, nTissue, domainID;} CHUNKS; 
typedef  struct column_st { int  nChunks; int offset ;} COLUMN; 
typedef  struct balance_domain_st { int id, nTissue, bbVol, bbHeight; } BALANCE_DOMAIN; 
typedef  struct balancer_st 
{ 
   int rank,  nRank, nTasks, nC, nx, ny, nz, dx, dy, dz, NX, NY, NXYz, nCellLocal; 
   long long unsigned nCellGlobal;
   double alphaWork; 
   double alphaBalance; 
   double beta; 
   double timeRef; 
   COLUMN *columns;
   CHUNKS *chunkBuffer;
   CHUNKS chunk; 
   unsigned short *segLocal; 
   unsigned short *seg; 
   int printStats; 
   MPI_Comm comm;
} BALANCER; 

#ifdef __cplusplus
extern "C" {
#endif
PARTITION_INFO  corePartition(CHUNKS chunk, BALANCER *balancer);
BALANCER buildBalancer(long long unsigned nCellGlobal, int nx,int ny,int nz,int dx,int dy,int dz, int  nTasks, int nCores, int nRCoresBB, double alpha, double beta, int printStat, MPI_Comm comm) ;
void fillSeg(int nn, long long unsigned *gidArray, BALANCER *balancer) ;
void reduceSeg(BALANCER *balancer) ;
int buildColumns(BALANCER *balancer);
BALANCE_DOMAIN domainInfo(long long unsigned  gid, BALANCER *balancer) ;
int balanceCores(BALANCER *balancer);
void  setDomainInfo(BALANCER *balancer);
void  printDomainInfo(BALANCER *balancer);
void destroyBalancer(BALANCER *balancer); 

#ifdef __cplusplus
}
#endif

#endif
