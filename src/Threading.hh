#ifndef THREADING_HH
#define THREADING_HH


class CoreGroup
{
 public:
   int threadID();
   int nThreads(){return nThreads_;}

   int nCores_;
   int firstCore_;
   int nThreads_;
   int* ompID_;
};

struct OMPTable { int omp_id, proc_id, core_id, hwThread_id;}; 
struct CoreTable { int core_id, first, nThread;}; 

void groupInfo(CoreGroup *group, int& coreID, int&hwThreadID, int& threadID, int& nCores, int& nHwThreads, int& nThreads) ;
class Threading 
{
 public:
   
   Threading(); 
   
   int nCores(){return nCores_;}
   int nThreads(){return nThreads_;}
   CoreGroup* mkGroup(int nCores) ; 
   
   CoreGroup** threadingMap_; 

 private:
   void probeHardware();
   
   CoreTable* coreTable_; 
   OMPTable*  ompTable_; 
   
   int nCores_;
   int nThreads_;
   int binding_; 
   int nRemainingCores_; 
   int nGroups_; 
   int nThreadsPerCore_;
   CoreGroup* groups_; 
};

#endif
