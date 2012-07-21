#ifndef WORKBOUND_H
#define WORKBOUND_H

typedef  struct chunks_st { int  zL, zU, nTissue, domainID;} CHUNKS; 
typedef  struct domain_st { int  nChunks; CHUNKS *chunks;} COLUMN; 
typedef  struct balancer_st { int rank,  nRank, nTasks, nx, ny, nz, dx, dy, dz, NX, NY, NXYz; CHUNKS *chunkBuffer;} BALANCER; 

#ifdef __cplusplus
extern "C" {
#endif

COLUMN *buildColumns(unsigned short *seg, BALANCER *balancer);
unsigned short *buildSeg(int nlocal, int *gidArray, BALANCER *balancer) ;
BALANCER buildGeom(int nx,int ny,int nz,int dx,int dy,int dz, int  nTasks) ;
int  domainID(int gid, COLUMN *column, BALANCER *balancer) ;

void destroyColumns(COLUMN *columns);
void destroySeg(unsigned short *seg);

#ifdef __cplusplus
}
#endif

#endif
