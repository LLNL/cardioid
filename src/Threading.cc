#include "Threading.hh"
#include <stdio.h> 
#include <stdlib.h> 
#include <omp.h> 
#include <cassert>
#include "mpiUtils.h"
#ifdef BGQ 
#include <spi/include/kernel/location.h> 
#else
#include <sys/types.h> 
#include <sys/sysctl.h> 
#endif
int cmpFuncOMPTABLE(const void *A, const void*B)
{
   OMPTABLE *tA = (OMPTABLE *)A;
   OMPTABLE *tB = (OMPTABLE *)B;
   if ( tA->core_id < tB->core_id)  return -1; 
   if ( tA->core_id > tB->core_id)  return  1; 
   if ( tA->hwThread_id < tB->hwThread_id)  return -1; 
   if ( tA->hwThread_id > tB->hwThread_id)  return  1; 
   return 0; 
}
void Threading::probeHardware()
{

   int omp_nThreads  = omp_get_max_threads(); 
   int omp_nProc = omp_get_num_procs(); 
   ompTable_=(OMPTABLE *)malloc(omp_nThreads*sizeof(OMPTABLE)); 
   #pragma omp parallel
   {
      int omp_id = omp_get_thread_num(); 
      int ncpu;
      #pragma omp critical
      {
         int proc_id, core_id, hwThread_id; 
         int flag =0; 
         #ifdef BGQ 
         proc_id     = Kernel_ProcessorID();       // 0-63
         core_id     = Kernel_ProcessorCoreID();   // 0-15
         hwThread_id = Kernel_ProcessorThreadID(); // 0-3
         flag =1; 
         #endif
         #ifdef __APPLE__
	 int requestNCPU[] = {CTL_HW,HW_NCPU}; 
	 int requestAVAILCPU[] = {CTL_HW,HW_AVAILCPU}; 
         ncpu,physicalCores,availcpu; 
         size_t len = sizeof(int); 
         sysctl(requestNCPU,2,&ncpu,&len,NULL,0); 
         sysctl(requestAVAILCPU,2,&availcpu,&len,NULL,0); 
         sysctlbyname("hw.physicalcpu", &physicalCores, &len, NULL, 0);

         proc_id = omp_id; 
         core_id = proc_id%ncpu;
         hwThread_id = 0; 
         core_id = omp_id%omp_nProc; 
         flag =1; 
         #endif
         if (!flag)
         {
            ncpu=omp_nProc; 
            proc_id=omp_id; 
            core_id = proc_id%ncpu;
            hwThread_id = 0; 
	 }
         ompTable_[omp_id].omp_id=omp_id; 
         ompTable_[omp_id].proc_id=proc_id; 
         ompTable_[omp_id].core_id=core_id; 
         ompTable_[omp_id].hwThread_id=hwThread_id; 
         int i = omp_id; 
      }
      #pragma omp barrier
      #pragma omp master
      {
        qsort(ompTable_,omp_nThreads,sizeof(OMPTABLE),cmpFuncOMPTABLE); 
        int core_id = -1; 
        int nCore=0; 
	for (int i =0;i<omp_nThreads;i++) if (ompTable_[i].core_id != core_id) nCore++; 
        coreTable_ = (CORETABLE *)malloc(nCore*sizeof(CORETABLE)); 
        nCore=0; 
        int first=0; 
	int last =0; 
        core_id = ompTable_[0].core_id; 
	for (int i =0;i<omp_nThreads;i++) 
        {
           if (ompTable_[i].core_id != core_id) 
           {
		printf("%d %d %d\n",core_id,first,last); 
                coreTable_[nCore].core_id=core_id; 
                coreTable_[nCore].first=first; 
                coreTable_[nCore].nThread=last-first+1; 
                first=i; core_id=ompTable_[i].core_id; 
                nCore++; 
           }
           last = i; 
        }
        coreTable_[nCore].core_id=core_id; 
        coreTable_[nCore].first=first; 
        coreTable_[nCore].nThread=last-first+1; 
        printf("%d %d %d\n",core_id,first,last); 
        nCore++; 
        nCores_=nCore; 
      }
      #pragma omp barrier
   }
    
}
Threading::Threading()
{
nThreads_  = omp_get_max_threads(); 
   probeHardware(); 
   nGroups_=0; 
   nRemainingCores_=nCores_; 
   groups_ = (coreGroup*)malloc(nCores_*sizeof(coreGroup)); 
   threadingMap_ =(coreGroup**)malloc(nThreads_*sizeof(coreGroup*)); 
   for (int ii=0;ii<nThreads_;ii++) threadingMap_[ii]=(coreGroup*)NULL; 
   
}
void groupInfo(coreGroup *group, int& coreID, int&hwThreadID, int& threadID, int& nCores, int& nHwThreads, int& nThreads) 
{
   threadID=groupThreadID(group);
   nThreads = group->nThreads;
   nCores= group->nCores; 
   nHwThreads=nThreads/nCores; 
   coreID = threadID/nHwThreads; 
   hwThreadID = threadID%nHwThreads; 
}

int groupThreadID(coreGroup *group) 
{
  int ompID = omp_get_thread_num(); 
  int threadID=-1; 
  for (int ii = 0;ii<group->nThreads;ii++) 
  {
     if ( ompID == group->ompID[ii] ) {threadID=ii; break;}
  }
  assert(threadID != -1) ; 
  return threadID; 
}
int Threading::nCores() { return nCores_ ; } 
coreGroup* Threading::mkGroup(int nCores ) 
{
   coreGroup group; 
   if (nCores==-1) nCores = nRemainingCores_; 
   group.nCores=nCores; 
   int firstCoreIndex=nCores_-nRemainingCores_; 
   int nThreads=0; 
   for(int ii=0;ii<nCores;ii++) nThreads+=coreTable_[ii+firstCoreIndex].nThread; 
   group.nThreads=nThreads;
   nRemainingCores_-=nCores; 
   group.ompID = (int *)malloc(group.nThreads*sizeof(int)); 
   assert(nRemainingCores_ >= 0); 
   int kk=0; 
   for (int ii=0;ii<nCores;ii++) 
   {
      int nThreads=coreTable_[ii+firstCoreIndex].nThread; 
      int first=coreTable_[ii+firstCoreIndex].first; 
      for(int jj=0;jj<nThreads;jj++) 
      {
         int ompID = ompTable_[first+jj].omp_id; 
       	 group.ompID[kk++]= ompID; 
       	 threadingMap_[ompID] = groups_+nGroups_; 
      }
   }
   groups_[nGroups_++] = group; 
   if (getRank(0)==0) printf("group=%p nCores=%d nThreads=%d",groups_+nGroups_-1,group.nCores,group.nThreads); 
   for(int ii=0;ii<group.nThreads;ii++) printf(" %d",group.ompID[ii]); printf("\n") ;
   return groups_+nGroups_-1;
}
#if (0) 
main()
{
   Threading info; 
   int nCores=info.nCores_; 
   coreGroup *groupA=info.mkGroup(1); 
   coreGroup *groupB=info.mkGroup(-1); 
}

#endif
