#include "ThreadServer.hh"
#include <algorithm>
#include <cassert>
#include <set>
#include <iostream>
#include <iomanip>
#ifdef BGQ 
#include <spi/include/kernel/location.h> 
#endif

using namespace std;


class CoreMatchPredicate
{
 public:
   CoreMatchPredicate(int coreID) : targetCore_(coreID){};
   bool operator()(const ThreadHardwareInfo& a){return targetCore_ == a.coreID_;}
 private:
   int targetCore_;
};

class CoreMatcher
{
 public:
   bool operator()(const ThreadHardwareInfo& a, const int& b)
   {return a.coreID_ < b;}
   bool operator()(const int& b, const ThreadHardwareInfo& a)
   {return b < a.coreID_;}
};

/** Initializes the hardware info for the calling thread.  This function
 *  is designed to be called by any and all threads */
ThreadHardwareInfo::ThreadHardwareInfo()
{
   // generic, non-harware specific
   int ncpu = omp_get_num_procs(); 
   //ncpu = 16; 
   ompID_ = omp_get_thread_num();
   procID_ = ompID_; 
   coreID_ = ompID_%ncpu;
   hwThreadID_ = ompID_/ncpu; 

   // The rest of this function is overrides for the generic info on the
   // specific hardware where better information is available.

   // BGQ 
   #ifdef BGQ
   procID_     = Kernel_ProcessorID();       // 0-63
   coreID_     = Kernel_ProcessorCoreID();   // 0-15
   hwThreadID_ = Kernel_ProcessorThreadID(); // 0-3
   #endif
   
   // Mac
   // Jim wrote some code that made sysctl calls to get information on
   // Mac, but none of that information was actually used.  For now,
   // there is no special Mac version.  To see what Jim did you'll need
   // to dig Threading.cc out of the repo from about r525.
}

void ThreadTeam::addMember(const ThreadHardwareInfo& thread)
{
   threads_.push_back(thread);
   buildRankInfo();
} 

void ThreadTeam::buildRankInfo()
{
   typedef multiset<ThreadHardwareInfo>::iterator Iter;
   typedef pair<Iter, Iter> IterPair;

   rankInfo_.resize(omp_get_max_threads());
   multiset<ThreadHardwareInfo> mset(threads_.begin(), threads_.end());

   nSquads_ = 0;
   int teamRank = 0;
   while (mset.size() > 0)
   {
      Iter iter = mset.begin();
      int targetCore = iter->coreID_;
      IterPair range = equal_range(mset.begin(), mset.end(), targetCore, CoreMatcher());
      int squadSize = distance(range.first, range.second);
      int squadRank = 0;
      for (Iter iter=range.first; iter!=range.second; ++iter)
      {
         rankInfo_[iter->ompID_].teamRank_  = teamRank++;
         rankInfo_[iter->ompID_].squadRank_ = squadRank++;
         rankInfo_[iter->ompID_].squadSize_ = squadSize;
         rankInfo_[iter->ompID_].coreRank_  = nSquads_;
      }
      ++nSquads_;
      mset.erase(range.first, range.second);
   }
}

ostream& operator<<(ostream& out, const ThreadTeam& tt)
{
   out << "nSquads=" << tt.nSquads() << " nThreads=" << tt.nThreads() << endl;
   out << "     Core hwThread   omp_id teamRank coreRank squadRank\n"
       << "--------------------------------------------------------" << endl;
   for (unsigned ii=0; ii<tt.nThreads(); ++ii)
   {
      unsigned ompID = tt.hwInfo(ii).ompID_;
      out << setw(9) << tt.hwInfo(ii).coreID_
          << setw(9) << tt.hwInfo(ii).hwThreadID_
          << setw(9) << tt.hwInfo(ii).ompID_
          << setw(9) << tt.rankInfo(ompID).teamRank_
          << setw(9) << tt.rankInfo(ompID).coreRank_
          << setw(9) << tt.rankInfo(ompID).squadRank_
          << endl;
   }
   return out;
}



bool operator<(const ThreadHardwareInfo& a, const ThreadHardwareInfo& b)
{
   if ( a.coreID_ < b.coreID_)  return true; 
   if ( a.coreID_ == b.coreID_ &&
        a.hwThreadID_ < b.hwThreadID_)  return true;
   return false;
}

ThreadServer& ThreadServer::getInstance()
{
   static ThreadServer instance;
   return instance;
}

/** Populate the list of available threads.  At the start, all threads
 * are available */
ThreadServer::ThreadServer()
{
   availableThreads_.reserve(omp_get_max_threads());
   #pragma omp parallel
   {
      // critical section since stl operations are not thread safe
      #pragma omp critical
      {
         availableThreads_.push_back(ThreadHardwareInfo());
      }
      #pragma omp barrier
   }

   sort(availableThreads_.begin(), availableThreads_.end());
}



ThreadTeam ThreadServer::getThreadTeam(const vector<unsigned>& coreID)
{
   ThreadTeam team;

   // special case: empty list.  Add all remaining cores to team.
   if (coreID.size() == 0)
   {
      // can't do dev with 1 thread per node.
      assert(availableThreads_.size() > 0); 
      for (unsigned ii=0; ii<availableThreads_.size(); ++ii)
         team.addMember(availableThreads_[ii]);
      availableThreads_.clear();
   }
   
   
   for (unsigned ii=0; ii<coreID.size(); ++ii)
   {
      CoreMatchPredicate coreMatch(coreID[ii]);
      vector<ThreadHardwareInfo>::iterator iter =
         find_if(availableThreads_.begin(), availableThreads_.end(), coreMatch);
      assert(iter != availableThreads_.end());
      team.addMember(*iter);
      availableThreads_.erase(iter);
   }

   return team;
}

