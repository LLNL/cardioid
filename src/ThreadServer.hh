#ifndef THREAD_SERVER_HH
#define THREAD_SERVER_HH

#include <vector>
#include <iosfwd>
#include <omp.h>

struct ThreadRankInfo
{
   ThreadRankInfo()
   : teamRank_(-1), squadRank_(-1), squadSize_(-1), coreRank_(-1){}
   int teamRank_;  // rank of thread in team
   int squadRank_; // rank of thread in squad
   int squadSize_;
   int coreRank_;  // rank of the core for this squad
};


class ThreadHardwareInfo
{
 public:
   ThreadHardwareInfo();
   int ompID_;
   int procID_;
   int coreID_;
   int hwThreadID_;
};

class ThreadTeam
{
 public:
   void addMember(const ThreadHardwareInfo& thread);
   int nThreads()  const {return threads_.size();}
   int nSquads()   const {return nSquads_;}
   int teamRank()  const {return rankInfo_[omp_get_thread_num()].teamRank_;}
   int teamRank(int ompId)  const {return rankInfo_[ompId].teamRank_;}
   int squadRank() const {return rankInfo_[omp_get_thread_num()].squadRank_;}
   int squadSize() const {return rankInfo_[omp_get_thread_num()].squadSize_;}
   const ThreadRankInfo& rankInfo()          const {return rankInfo_[omp_get_thread_num()];}
   const ThreadRankInfo& rankInfo(int ompID) const {return rankInfo_[ompID];}
   const ThreadHardwareInfo& hwInfo(int ii) const {return threads_[ii];} // not ranked by ompId or teamRank
   ThreadTeam() { nSquads_ = 0; }
   
 private:
   void buildRankInfo();
   int nSquads_;
   std::vector<ThreadHardwareInfo> threads_; // info only on threads in team
   std::vector<ThreadRankInfo> rankInfo_;    // info for all omp ids.
};

std::ostream& operator<<(std::ostream& out, const ThreadTeam& tt);


/** This implementation isn't even remotely thread safe.  Only call from
 *  unthreaded regions of the code.
 *
 *  This class serves requests for thread teams by keeping track of
 *  which threads have already been requested and ensuring that no
 *  thread is assigned to multiple teams
 */
class ThreadServer
{
 public:
   static ThreadServer& getInstance();
   ThreadTeam getThreadTeam(const std::vector<unsigned>& coreID);
   
 private:
   ThreadServer();
   ThreadServer(const ThreadServer&);
   ThreadServer& operator=(const ThreadServer&);
   
   std::vector<ThreadHardwareInfo> availableThreads_;
};

#endif
