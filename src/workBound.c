#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

typedef  struct chunks_st { int  zL, zU, nTissue, domainID;} CHUNKS; 
typedef  struct domain_st { int  nChunks; CHUNKS *chunks;} COLUMN; 
typedef  struct balancer_st { int rank,  nRank, nTasks, nx, ny, nz, dx, dy, dz, NX, NY, NXYz; CHUNKS *chunkBuffer;} BALANCER; 

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
  int task=0; 
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
    if (chunkBuffer != NULL) chunks = chunkBuffer+task;  else chunks=NULL; 
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
              chunks[nChunks].domainID=task; 
              }
              nChunks++; 
              task++; 
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
          chunks[nChunks].domainID=task; 
       }
       nChunks++; 
       task++; 
     }
     if (column != NULL) 
     {
       column[cxy].nChunks = nChunks; 
       if (nChunks > 0) column[cxy].chunks =  chunks; 
       else column[cxy].chunks=NULL; 
       sumChunks += nChunks; 
     }
   }
   printf("sumChunks=%d\n",sumChunks); 
   return task; 
}
void allocChunkBuffer(BALANCER *balancer)
{
   balancer->chunkBuffer = (CHUNKS *)malloc(sizeof(CHUNKS)*balancer->nTasks); 
}
void destroyColumns(COLUMN *columns) 
{ 
   free(columns); 
}
COLUMN *buildColumns(unsigned short *seg,BALANCER *balancer)
{
   int NX = balancer->NX; 
   int NY = balancer->NY; 
   COLUMN *columns = (COLUMN *)malloc(sizeof(COLUMN)*NX*NY); 
   allocChunkBuffer(balancer); 
   if (balancer->rank==0) 
   {
      int nx = balancer->nx; 
      int ny = balancer->ny; 
      int nz = balancer->nz; 
      int dx = balancer->dx; 
      int dy = balancer->dy; 
      int dz = balancer->dz; 
      double minA = -1; 
      double maxA = -1; 
      double  a ; 
      int nTasks; 
      for (a = 0.0125;a<2;a*=2.0) 
      {
         nTasks=findChunks(balancer,seg, a,NULL);
         if ( nTasks < balancer->nTasks) minA = a; 
         if ( nTasks > balancer->nTasks) {maxA = a; break;}
         printf("min/max %f %f %f %d %d\n",a,minA,maxA,nTasks,balancer->nTasks); 
      }
      printf("min/max %f %f %f %d %d\n",a,minA,maxA,nTasks,balancer->nTasks); 
      int loop=0; 
      while(nTasks != balancer->nTasks   )
      {
        printf("min/max %f %f %f %d\n",a,minA,maxA,nTasks); 
        a = 0.5*(minA+maxA); 
        nTasks=findChunks(balancer,seg, a, NULL);
        if (nTasks < balancer->nTasks) minA = a; 
        if (nTasks > balancer->nTasks) maxA = a; 
        if (loop++ > 30) break; 
      }
      printf("min/max %f %f %f %d\n",a,minA,maxA,nTasks); 
      if (nTasks == balancer->nTasks)  
      {
         nTasks=findChunks(balancer,seg, a, columns);
      }
   }
   MPI_Bcast(columns,NX*NY*sizeof(COLUMN),MPI_BYTE,0,MPI_COMM_WORLD); 
   MPI_Bcast(balancer->chunkBuffer,balancer->nTasks*sizeof(CHUNKS),MPI_BYTE,0,MPI_COMM_WORLD); 

   CHUNKS *chunkBufferRoot = balancer->chunkBuffer; 
   MPI_Bcast(&(chunkBufferRoot),sizeof(CHUNKS *),MPI_BYTE,0,MPI_COMM_WORLD); 
   assert(chunkBufferRoot!=NULL); 
   long long correction = balancer->chunkBuffer-chunkBufferRoot; 

   for (int l=0;l<NX*NY;l++) if (columns[l].nChunks > 0) columns[l].chunks += correction; 
     
   MPI_Barrier(MPI_COMM_WORLD); 
   return columns;
}
int  domainID(int gid, COLUMN *column, BALANCER *balancer) 
{
   int x=gid%balancer->nx; 
   int y=((gid-x)/balancer->nx)%balancer->ny; 
   int z=((gid-x-balancer->nx*y))/(balancer->nx*balancer->ny); 
   int i = x/balancer->dx;
   int j = y/balancer->dy;
   int cxy = (i + balancer->NX*j); 
   for (int ic=0;ic<column[cxy].nChunks;ic++) 
      if (column[cxy].chunks[ic].zL <= z &&  z<= column[cxy].chunks[ic].zU) return column[cxy].chunks[ic].domainID; 
   assert(0); 
   return -1; 
}
void destroySeg(unsigned short *seg) { free(seg);}
unsigned short *buildSeg(int nlocal, int *gidArray, BALANCER *balancer) 
{
   int NX = balancer->NX; 
   int NY = balancer->NY; 
   int NXYz =balancer->NXYz;
   int nx = balancer->nx; 
   int ny = balancer->ny; 
   int nz = balancer->nz; 
   int dx = balancer->dx; 
   int dy = balancer->dy; 
   int dz = balancer->dz; 
   unsigned short *seg = (unsigned short *)malloc(sizeof(*seg)*NXYz); 
   for (int i=0;i<NXYz;i++) seg[i]=0; 
   int minX,minY,minZ;
   int minXG,minYG,minZG;
   int maxX,maxY,maxZ;
   int maxXG,maxYG,maxZG;
   for (int ii=0;ii<nlocal;ii++) 
   {
      unsigned  gid = gidArray[ii]; 
      int x=gid%nx; 
      int y=((gid-x)/nx)%ny; 
      int z=((gid-x-nx*y))/(nx*ny); 
      int i = x/dx;
      int j = y/dy;
      int segIndex = z + nz*(i + NX*j); 
      seg[segIndex] += 1; 
      if (i==0 || x < minX) minX = x; 
      if (i==0 || y < minY) minY = y; 
      if (i==0 || z < minZ) minZ = z; 
      if (i==0 || x > maxX) maxX = x; 
      if (i==0 || y > maxY) maxY = y; 
      if (i==0 || z > maxZ) maxZ = z; 
   }
   MPI_Allreduce(&minX,&minXG,1,sizeof(int),MPI_MIN,MPI_COMM_WORLD); 
   MPI_Allreduce(&minY,&minYG,1,sizeof(int),MPI_MIN,MPI_COMM_WORLD); 
   MPI_Allreduce(&minZ,&minZG,1,sizeof(int),MPI_MIN,MPI_COMM_WORLD); 
   MPI_Allreduce(&maxX,&maxXG,1,sizeof(int),MPI_MAX,MPI_COMM_WORLD); 
   MPI_Allreduce(&maxY,&maxYG,1,sizeof(int),MPI_MAX,MPI_COMM_WORLD); 
   MPI_Allreduce(&maxZ,&maxZG,1,sizeof(int),MPI_MAX,MPI_COMM_WORLD); 
   if (balancer->rank==0) 
   {
      printf("min %d %d %d\n",minXG,minYG,minZG); 
      printf("max %d %d %d\n",maxXG,maxYG,maxZG); 
   }
   MPI_Barrier(MPI_COMM_WORLD); 
   
   unsigned short *segtmp = (unsigned short *)malloc(sizeof(*seg)*NXYz); 
   for (int i=0;i<NXYz;i++) segtmp[i] = seg[i]; 
   MPI_Reduce(segtmp,seg,NXYz,MPI_UNSIGNED_SHORT,MPI_SUM,0,MPI_COMM_WORLD); 
   free(segtmp); 
   return seg; 
}
BALANCER buildGeom(int nx,int ny,int nz,int dx,int dy,int dz, int  nTasks) 
{
   int rank,nRank; 
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nRank);
   BALANCER balancer; 
   balancer.rank = rank; 
   balancer.nRank = nRank; 
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
   balancer.chunkBuffer=NULL; 
   return balancer; 
}
#if 0 
int main(int argc, char *argv[])
{
   MPI_Init(&argc,&argv); 
   int rank,nRank; 
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nRank);
   int nx=960; 
   int ny=1050; 
   int nz=930; 
   int dx=18;
   int dy=16;
   int dz=14; 
   unsigned int *buff = (unsigned int *)malloc(sizeof(unsigned int)*380000000); 
   char filename[32]; 
   int ngrid=0; 
   int n; 
   for (int ii=rank;ii<rank+4/nRank;ii++) 
   {
   sprintf(filename,"gid#%6.6d",ii); 
   printf("%s\n",filename); 
   FILE *file = fopen(filename,"r"); 
   size_t size = 16*1024 ;
   while(!feof(file))
   {
       n = fread(buff+ngrid,sizeof(*buff),size,file);   // In Cardiod just look over cells on the node. 
       ngrid+=n;
   }
   fclose(file); 
   printf("n=%d %d\n",n,ngrid); 
   }

     
   int target=96*1024; 

   BALANCER balancer =buildGeom(nx,ny,nz,dx,dy,dz,target); 
   unsigned short *seg=buildSeg(ngrid,buff,&balancer); 
   COLUMN *columns = buildColumns(seg,&balancer);
   destroySeg(seg); 
   sprintf(filename,"Stuff#%6.6d",balancer.rank); 
   FILE *stuff = fopen(filename,"w"); 
   for (int i =0;i<ngrid;i++) 
   {
    int gid = buff[i]; 
    int id = domainID( gid, columns, &balancer) ;
    fprintf(stuff,"%d %d\n",gid,id); 
   }
   if (rank == 0) destroyColumns(columns); 
}
#endif
