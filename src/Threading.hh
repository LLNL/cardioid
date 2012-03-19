#ifndef THREADING_HH
#define THREADING_HH

#include <vector>

/** I would prefer that this class not have a default constructor.
 * However, I'm forced into if right now by the requirement that the
 * Threading class has groups_ as a c-array of CoreGroup.  If we call
 * malloc to get memory for that array bad things happen since the
 * memory will not be initialized and trying to assign into the array
 * will cause havoc.  To use new for the array requires a default
 * constructor.
 *
 *  Eventually I'll figure out how to get the pointers out of the
 *  threading class.  Then I can get rd of the default constructor.
 */
class CoreGroup
{
 public:
   CoreGroup(int nCores, int nThreads, const std::vector<int>& ompID)
   : nCores_(nCores), nThreads_(nThreads), ompID_(ompID){}
   CoreGroup()
   : nCores_(0), nThreads_(0), ompID_(std::vector<int>(0)){}
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
