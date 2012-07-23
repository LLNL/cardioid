#ifndef WORKBOUND_H
#define WORKBOUND_H
#include <mpi.h>

typedef  struct chunks_st { int  zL, zU, nTissue;} CHUNKS; 
typedef  struct column_st { int  nChunks; int offset ;} COLUMN; 
typedef  struct balance_domain_st { int id, nTissue, bbVol, bbHeight; } BALANCE_DOMAIN; 
typedef  struct balancer_st { int rank,  nRank, nTasks, nx, ny, nz, dx, dy, dz, NX, NY, NXYz; 
                                 double relCost; 
                                 COLUMN *columns;
                                 CHUNKS *chunkBuffer;
                                 unsigned short *segLocal; 
                                 unsigned short *seg; 
                                 int printStats; 
                                 MPI_Comm comm;
                            } BALANCER; 

#ifdef __cplusplus
extern "C" {
#endif
BALANCER buildBalancer(int nx,int ny,int nz,int dx,int dy,int dz, int  nTasks, int printStat, MPI_Comm comm) ;
void fillSeg(int nn, long long unsigned *gidArray, BALANCER *balancer) ;
void reduceSeg(BALANCER *balancer) ;
int buildColumns(BALANCER *balancer);
BALANCE_DOMAIN domainInfo(long long unsigned  gid, BALANCER *balancer) ;
void  printDomainInfo(BALANCER *balancer);
void destroyBalancer(BALANCER *balancer); 

#ifdef __cplusplus
}
#endif

#endif
