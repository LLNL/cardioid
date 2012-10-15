#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <unistd.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h> 
#include <string.h> 
#include "TT06Gates.h"
void HPM_Start(char *); 
void HPM_Stop(char *); 
typedef struct { int nCores, nThreads; struct { int coreID, threadID;} map[64];} threadInfo; 
//#define BGQ
#ifdef BGQ
#include <spi/include/l1p/pprefetch.h> 
#include <spi/include/l1p/sprefetch.h> 
#endif 
#include "TT06Func.h"
void update_nonGate(void *ptr, double dt, struct CellTypeParms *p, int nCells, int *cellTypeVector, const double *VM, int offset, double *gates[19], double *dVdt);
void update_nonGate_v1(void *ptr, double dt, struct CellTypeParms *p,int nCells, int *cellTypeVector, const double *VM, int offset, double *gates[19], double *dVdt);
double initState(double *states,double *gates, int cellType);
char *getStateName(int index);


void initExp(); 
void fastLogInit();
//void initArray();

typedef void (*UPDATENONGATE)(void *ptr, double dt, struct CellTypeParms *p,int nCells, int *cellTypeVector, const double *VM, int offset, double *gates[19], double *dVdt);
typedef void  (*UPDATEGATE)(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a) ; 

static   double mhu[13][50]; 
static   double tauR[13][50]; 
static   double dt=0.01; 
static   int nonGatesFlag=0;
static   int gatesFlag=0;

threadInfo getThreadInfo()
{
   int nThreads = omp_get_max_threads(); 

   int nCores= 4; 
#ifdef BGQ
    nCores=16; 
#endif 
   threadInfo info ;
   info.nCores=nCores; 
   info.nThreads=nThreads; 
   for (int ompID=0;ompID<nThreads;ompID++) 
   {
	info.map[ompID].coreID = ompID%nCores; 
	info.map[ompID].threadID =  (ompID/nCores) *nCores  + info.map[ompID].coreID;
   }
   return info; 
}
   
void init(int nCellsOnNode, int *cellTypeVector, double **g, double *Vm)
{
      int ompID = omp_get_thread_num(); 
   for (int i=0;i<nCellsOnNode;i++) 
   {
        
        double state[nStateVar]; 
        double gate[1]; 
        int cellType=0; 
        double vm=initState(state,gate,cellType); 
        cellTypeVector[i] = cellType; 
        Vm[i] = vm; 
	     for (int j=0;j<19;j++)  g[j][i]=state[j];
	     g[19][i]=state[18];
   }
}
void readData()
{
    
   char line[1024],name[16]; 
   int l,m; 
   FILE *file=fopen("../coef.data","r"); 
   for (int i=0;i<13;i++)
   {
      double *a; 
      int k; 
      fgets(line,1023,file); 
      a = mhu[i]; 
      k=sscanf(line,"%s %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf",
      name,&l,&m,a+0,       a+1,a+2,a+3,a+4,a+5,a+6,a+7,a+8,a+9,a+10,a+11,a+12,a+13,a+14,a+15,a+16,a+17,a+18,a+19,a+20,a+21,a+22,a+23,a+24); 

      fgets(line,1023,file); 
      a = tauR[i]; 
      k=sscanf(line,"%s %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf",
      name,&l,&m,a+0,       a+1,a+2,a+3,a+4,a+5,a+6,a+7,a+8,a+9,a+10,a+11,a+12,a+13,a+14,a+15,a+16,a+17,a+18,a+19,a+20,a+21,a+22,a+23,a+24); 
   }
   fclose(file); 
}
void parallelSection(threadInfo *info,  UPDATENONGATE updateNonGateFunc, UPDATEGATE updateGateFuncs[], double dt, int maxLoop, int nCells, int *cellTypeVector, double *Vm, double *g[19], double *dVdt, FILE *file)
{
  struct CellTypeParms cellParms[4]={{0,0,2.724,0.392,0.073,0.0},{0,0,2.724,0.098,0.294,0.0},{ 0,0,2.724, 0.392, 0.294, 0.0 },{0,0,0.2724,0.098,0.294,0.0}}; 
#ifdef BGQ
   HPM_Start("Reaction"); 
   vprof_start(); 
#endif 
   double cpuTimes[128]; 
   for(int i=0;i<128;i++) cpuTimes[i]=0.0; 
   int loop = 0; 
   
   #pragma omp parallel
   {

      int nThreads   = info ->nThreads; 
      int nCores     = info ->nCores; 
      int nSquad = nThreads/nCores; 

      int nCellPerCore = nCells/nCores; 
      int nCellPerThread = nCells/nThreads; 

      int ompID = omp_get_thread_num(); 

      int coreID     = info->map[ompID].coreID; 
      int offsetCore =  coreID*nCellPerCore; 

      int threadID   = info->map[ompID].threadID; 
      int offsetThread =  threadID*(nCellPerThread); 

      int squadID = threadID / nCores; 
      int nEq = 12/nSquad;
      int offsetGEq =   squadID* nEq;
      
      int squadRank = threadID % nSquad; 

      int loop;  
      double time=0.0; 
      int map[12]={0,1,2,3,4,5,6,7,8,9,10,11}; 
      for (loop=0;loop<maxLoop;loop++)
      {
             double t0,t1; 
             if (nonGatesFlag)
             {
             #pragma omp barrier 
             t0 = omp_get_wtime(); 
             updateNonGateFunc(NULL,dt, cellParms, nCellPerThread, cellTypeVector, Vm+offsetThread, offsetThread, g, dVdt+offsetThread);
             t1 = omp_get_wtime(); 
             cpuTimes[2*ompID]+= t1-t0; 
             }

             if (gatesFlag)
             {
             #pragma omp barrier 
             t0 = omp_get_wtime(); 
             for (int i=0;i<nEq;i++) 
             {
               int eqx = offsetGEq+i; 
               int eq = map[eqx]; 
                 updateGateFuncs[eq](dt, nCellPerCore, Vm+offsetCore, g[eq+7]+offsetCore, mhu[eq], tauR[eq]);  
             }
             t1 = omp_get_wtime(); 
             cpuTimes[2*ompID+1]+= t1-t0; 
             }
             time+=dt; 
      }
     #pragma omp barrier 
   }
#ifdef BGQ
   vprof_stop(); 
   HPM_Stop("Reaction"); 
#endif
   int nThreads   = info ->nThreads; 
   for (int i=0;i<19;i++) printf("%4d %4d %3d %24.14e %24.14e %24.14e\n", maxLoop,nCells,i,g[i][0],g[i][nCells/2],g[i][nCells-1] );
   if (file != NULL) 
   {
   fprintf(file,"#ompID coreID nonGate    gate\n");
   for(int ompID=0;ompID<nThreads;ompID++) fprintf(file,"%6d %6d %8.3f %8.3f\n",ompID,ompID%16,cpuTimes[2*ompID],cpuTimes[2*ompID+1]);fflush(stdout);
   fprintf(file,"end_of_data\n");
   }


}
int main(int argc, char **argv)
{
   int maxLoop = 100000; 
   int nCellsOnNode=4096;
   int jobId=Kernel_GetJobID(); 
   printf("jobId=%d\n",jobId); 
   int mode = 0775;
   char dirname[256]; 
   sprintf(dirname,"jobId.%d",jobId);
   int rc; 
   rc = mkdir(dirname, mode);
   rc = chdir(dirname); 

   FILE *file=fopen("time.data","w"); 
   for (int i=1;i<argc;i++)
   {
      if (strcmp(argv[i],"-nonGates") ==0)  nonGatesFlag=1; 
      if (strcmp(argv[i],"-gates") ==0)  gatesFlag=1; 
      if (strcmp(argv[i],"-n") ==0)  nCellsOnNode = atol(argv[++i]); 
      if (strcmp(argv[i],"-maxLoop") ==0)  maxLoop = atol(argv[++i]); 
      if (strcmp(argv[i],"-h") ==0)  
      {
	printf("usage: driver [-nonGates] [-gates] [-maxLoop <#TimeSteps>] [-n <#Particles>]\n"); 
        exit(0); 
      }
   }
   printf("%d %d %d %d\n",nonGatesFlag,gatesFlag,nCellsOnNode,maxLoop);  
   int myid; 
   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid); 
   
   UPDATEGATE updateFuncs0[]={ update_mGate, update_hGate, update_jGate, update_Xr1Gate, update_Xr2Gate, update_XsGate, update_rGate, update_dGate, 
                          update_fGate, update_f2Gate, update_jLGate, update_s0Gate, update_s1Gate} ;
   UPDATEGATE updateFuncs1[]={ update_mGate_v1, update_hGate_v1, update_jGate_v1, update_Xr1Gate_v1, update_Xr2Gate_v1, update_XsGate_v1, update_rGate_v1, update_dGate_v1, 
                          update_fGate_v1, update_f2Gate_v1, update_jLGate_v1, update_s0Gate_v1, update_s1Gate_v1} ;

// create a aligned buffer and offset 
   uint64_t sizeBuffer   = 22*(nCellsOnNode+63)*sizeof(double); 
   uint64_t buffer  = (uint64_t)malloc(sizeBuffer); 
   uint64_t  offset = 4*((nCellsOnNode+3)/4); 
   double *start  = (double *)(32*((buffer+31)/32));
   
//  Make sure the gates (g), nonGate (cell),  voltage (Vm) and dVm/dt (dVdt) are 32 byte aligned 

   double *g[20]; for (int i=0;i<20;i++) g[i] = start + i*offset;
   double *Vm =start+20*offset; 
   double *dVdt =start+21*offset; 
   int cellTypeVector[offset]; 
#if BGQ
#pragma omp parallel
{
  L1P_SetStreamPolicy(L1P_stream_confirmed);
  L1P_SetStreamAdaptiveMode(0);
  L1P_SetStreamDepth(2);   
} 
#endif 
   initExp(); 
   fastLogInit(); 
   initCnst(); 
   initNonGateCnst(); 
   readData();

   threadInfo info = getThreadInfo(); 
   
   init(nCellsOnNode, cellTypeVector,g, Vm);

   parallelSection(&info, update_nonGate, updateFuncs0, dt, maxLoop, nCellsOnNode, cellTypeVector, Vm, g, dVdt,NULL);
   double sum0[19]; 
   for (int eq =0;eq<19;eq++)
   {
      sum0[eq]=0; 
      for (int i=0;i<nCellsOnNode;i++) sum0[eq]+=g[eq][i];  sum0[eq] /= nCellsOnNode; 
   }

   init(nCellsOnNode, cellTypeVector,g, Vm);
   parallelSection(&info, update_nonGate_v1, updateFuncs1, dt, maxLoop, nCellsOnNode, cellTypeVector, Vm, g, dVdt,file);

   printf("\n********************************\n"); 
   printf("g0 = gSNorm[eq][0]\n");    
   printf("g1 = gSimd[eq][0]\n");    
   printf("ave  = 1 -0.5*(<g0>+<g1>)/g1\n"); 
   printf("diff  = <g1>-<g1>)/<g0>\n"); 
   printf("%2s %-9s %9s %9s %9s %15s %9s\n","eq","eq Name","g1","<g0>","<g1>","ave  ","diff");    
   for (int eq =0;eq<19;eq++)
   {
     double sum1=0; for (int i=0;i<nCellsOnNode;i++)   sum1+=g[eq][i]; sum1 /= nCellsOnNode; 
     double sum = g[eq][0]; 
     printf("%2d %-9s %10.3e %10.3e %10.3e error=%9.2e %9.2e\n",eq,getStateName(eq), sum,sum0[eq],sum1,1.0-0.5*(sum0[eq]+sum1)/sum,(sum1-sum0[eq])/sum0[eq]);    
   }

   MPI_Finalize(); 
}
