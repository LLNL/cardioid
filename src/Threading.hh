#ifndef THREADING_HH
#define THREADING_HH
using namespace std;
enum threadBinding { ROUNDROBIN} ;
typedef struct { int nCores; int firstCore; int nThreads; int *ompID; } coreGroup;
typedef  struct ompTable_st { int omp_id, proc_id, core_id, hwThread_id;} OMPTABLE; 
typedef  struct coreTable_st { int core_id, first, nThread;} CORETABLE; 
int groupThreadID(coreGroup *group) ;
void groupInfo(coreGroup *group, int& coreID, int&hwThreadID, int& threadID, int& nCores, int& nHwThreads, int& nThreads) ;
class Threading 
{
 public:

   Threading(); 
   void probeHardware(); 
   CORETABLE *coreTable_; 
   OMPTABLE  *ompTable_; 

   int nCores_;
   int nThreads_;
   int binding_; 
   int nRemainingCores_; 
   int nGroups_; 
   int nThreadsPerCore_;
   coreGroup **threadingMap_; 
   coreGroup *groups_; 
   int nCores();
   coreGroup* mkGroup(int nCores) ; 
};

#endif
