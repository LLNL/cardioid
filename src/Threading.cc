#include "Threading.hh"
#include <stdio.h> 
#include <omp.h> 
#include <cassert>
static int _nHwThreads;
Threading::Threading()
{
   nCores_    = omp_get_num_procs(); 
   nThreads_  = omp_get_max_threads(); 
   binding_   = ROUNDROBIN; 
   if ( nThreads_ < nCores_)
   {
     nThreadsPerCore_=1; 
     nCores_=nThreads_; 
   }
   else
   {
      assert (nThreads_%nCores_ == 0) ; 
      nThreadsPerCore_= nThreads_/nCores_; 
   }
   _nHwThreads = nThreadsPerCore_; 
   printf("nCores=%d nThreads=%d\n",nCores_,nThreads_); fflush(stdout); 
   nGroups_=0; 
   nRemainingCores_=nCores_; 
   groups_ = (coreGroup*)malloc(nCores_*sizeof(coreGroup)); 
   threadingMap_ =(coreGroup**)malloc(nThreads_*sizeof(coreGroup*)); 
   for (int ii=0;ii<nThreads_;ii++) threadingMap_[ii]=(coreGroup*)NULL; 
   
}
void groupInfo(coreGroup *group, int& coreID, int&hwThreadID, int& threadID, int& nCores, int& nHwThreads, int& nThreads) 
{ 
  int ompID = omp_get_thread_num(); 
  nHwThreads  = omp_get_max_threads()/omp_get_num_procs(); 
  nHwThreads  = _nHwThreads; 
  threadID=-1; 
  coreID = -1; 
  nThreads=group->nThreads; 
  nCores = group->nCores; 
  for (int ii = 0;ii<group->nThreads;ii++) 
  {
     if ( ompID == group->ompID[ii] ) {threadID=ii; break;}
  }
  assert(threadID != -1) ; 
  coreID = threadID/nHwThreads; 
  hwThreadID = threadID%nHwThreads; 
} 
int groupThreadID(coreGroup *group) 
{
int coreID=0,hwThreadID=0,threadID=0,nCores=0,nHwThreads=0,nThreads=0; 
groupInfo(group, coreID, hwThreadID, threadID, nCores, nHwThreads, nThreads) ;
return threadID; 
}
int Threading::nCores() { return nCores_ ; } 
coreGroup* Threading::mkGroup(int nCores ) 
{
   coreGroup group; 
   if (nCores==-1) nCores = nRemainingCores_; 
   group.nCores=nCores; 
   group.nThreads=nThreadsPerCore_*nCores;
   group.firstCore=nCores_-nRemainingCores_; 
   nRemainingCores_-=nCores; 
   group.ompID = (int *)malloc(group.nThreads*sizeof(int)); 
   assert(nRemainingCores_ >= 0); 
   int kk=0; 
   for (int ii=0;ii<nCores;ii++) 
   {
      int core = ii + group.firstCore; 
      for(int jj=0;jj<nThreadsPerCore_;jj++) 
      {
         int ompID = core + jj*nCores_; 
       	 group.ompID[kk++]= ompID; 
       	 threadingMap_[ompID] = groups_+nGroups_; 
      }
   }
   groups_[nGroups_++] = group; 
   printf("group=%p nCores=%d nThreads=%d",groups_+nGroups_-1,group.nCores,group.nThreads); 
   for(int ii=0;ii<group.nThreads;ii++) printf(" %d",group.ompID[ii]); printf("\n") ;
   return groups_+nGroups_-1;
}
/*
main()
{
   Threading info; 
   int nCores=info.nCores_; 
   coreGroup *groupA=info.mkGroup(1); 
   coreGroup *groupB=info.mkGroup(-1); 
}
*/
