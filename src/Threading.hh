#ifndef THREADING_HH
#define THREADING_HH
using namespace std;
enum threadBinding { ROUNDROBIN} ;
typedef struct { int nCores; int firstCore; int nThreads; int *ompID; } coreGroup;
int groupThreadID(coreGroup *group) ;
void groupInfo(coreGroup *group, int& coreID, int&hwThreadID, int& threadID, int& nCores, int& nHwThreads, int& nThreads) ;
class Threading 
{
 public:

   Threading(); 
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
