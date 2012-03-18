#ifndef THREADING_HH
#define THREADING_HH

#include <vector>

class CoreGroup
{
 public:
   CoreGroup(int nCores, int nThreads, const std::vector<int>& ompID)
   : nCores_(nCores), nThreads_(nThreads), ompID_(ompID){}
   int threadID() const;
   int nThreads() const {return nThreads_;}
   int nCores()   const {return nCores_;}

   std::vector<int> ompID_;
 private:
   int nCores_;
//   int firstCore_;
   int nThreads_;
};


void groupInfo(CoreGroup *group, int& coreID, int&hwThreadID, int& threadID, int& nCores, int& nHwThreads, int& nThreads) ;

struct OMPTable { int omp_id, proc_id, core_id, hwThread_id;}; 

class Threading 
{
   struct CoreTable { int core_id, first, nThread;}; 

 public:
   
   Threading(); 
   
   int nCores(){return nCores_;}
   int nThreads(){return nThreads_;}
   CoreGroup* mkGroup(int nCores) ; 
   
   CoreGroup** threadingMap_; 

 private:
   void probeHardware();
   
   std::vector<CoreTable> coreTable_; 
   std::vector<OMPTable>  ompTable_; 
   
   int nCores_;
   int nThreads_;
   int binding_; 
   int nRemainingCores_; 
   int nGroups_; 
   int nThreadsPerCore_;
   CoreGroup* groups_; 
};

#endif
