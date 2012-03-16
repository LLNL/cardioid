#include "PerformanceTimers.hh"
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <sys/time.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <omp.h>
#include "pio.h"
#include "ioUtils.h"
#include "tagServer.h"
#include "unionOfStrings.hh"

using namespace std;
namespace PerformanceTimers
{
enum TimerEnum { NCALLS, CYCLES} ;
/*
#ifdef BGL
uint64_t getTime()
{
   unsigned long long rts_get_timebase(void);
   uint64_t  t = 0.0;
   t = rts_get_timebase();
   return t; 
}
double getTick()
{
   double seconds_per_cycle = 1.4285714285714285714e-9;
   return seconds_per_cycle; 
}
#endif
*/

#include "PerformanceTimersBGQ.hh"
#include "PerformanceTimersGeneric.hh"

   struct TimerStruct
   {
      uint64_t start[8];
      uint64_t total[8];
   };
   TimerHandle loopIOTimer;
   TimerHandle sensorTimer;
   TimerHandle haloTimer;
   TimerHandle diffusionTimer;
   TimerHandle stimulusTimer;
   TimerHandle reactionTimer;
   TimerHandle integratorTimer;
   TimerHandle nonGateTimer;
   TimerHandle gateTimer;
   TimerHandle diffusionLoopTimer;
   TimerHandle reactionLoopTimer;
   TimerHandle reactionWaitTimer;
   TimerHandle diffusionWaitTimer;
   TimerHandle diffusionStallTimer;
   TimerHandle diffusionImbalanceTimer;
   TimerHandle imbalanceTimer;
   TimerHandle dummyTimer;
   
   
   vector<TimerStruct> timers_;
   typedef map<string, TimerHandle> HandleMap;
   HandleMap handleMap_;
   vector<string> printOrder_;
   string refTimer_;

}

using namespace PerformanceTimers;

void  profileInit()
{
   loopIOTimer = profileGetHandle("LoopIO");
   sensorTimer = profileGetHandle("Sensors");
   haloTimer = profileGetHandle("Halo Exchange");
   diffusionTimer= profileGetHandle("Diffusion");
   stimulusTimer = profileGetHandle("Stimulus");
   reactionTimer= profileGetHandle("Reaction");
   nonGateTimer= profileGetHandle("Reaction_nonGate");
   gateTimer = profileGetHandle("Reaction_Gate");
   diffusionLoopTimer= profileGetHandle("DiffusionLoop");
   integratorTimer = profileGetHandle("Integrator");
   haloTimer = profileGetHandle("Halo Exchange");
   reactionLoopTimer = profileGetHandle("ReactionLoop");
   reactionWaitTimer = profileGetHandle("ReactionWait");
   diffusionWaitTimer = profileGetHandle("DiffusionWait");
   diffusionStallTimer = profileGetHandle("DiffusionStall");
   diffusionImbalanceTimer = profileGetHandle("DiffusionImbalance");
   imbalanceTimer = profileGetHandle("Imbalance");
   dummyTimer = profileGetHandle("Dummy");
   machineSpecficInit(); 
}
       
void profileStart(const TimerHandle& handle)
{
   int tid=omp_get_thread_num() ;
   int id=handle+tid; 
   timers_[id].start[CYCLES] = getTime();
   for (int i=0; i<nCounters_;i++) 
   {
      uint64_t counter;
      readCounter(counterHandle[tid], i, &counter);   
      timers_[id].start[i+1] = counter; 
   }
}

void profileStop(const TimerHandle& handle)
{
   int tid=omp_get_thread_num() ;
   int id=handle+tid; 
   timers_[id].total[NCALLS] += 1;
   uint64_t delta = getTime() - timers_[id].start[CYCLES];
   timers_[id].total[CYCLES] += delta;
   for (int i=0; i<nCounters_;i++) 
   {
      uint64_t counter;
      readCounter(counterHandle[tid], i, &counter);
      uint64_t delta  = counter - timers_[id].start[i+1];
      timers_[id].total[i+1] += delta;
   }
}

void profileStart(const std::string& timerName)
{
   profileStart(profileGetHandle(timerName));
}

void profileStop(const std::string& timerName)
{
   profileStop(profileGetHandle(timerName));
}

TimerHandle profileGetHandle(string timerName)
{
    TimerHandle handle;
    TimerHandle handleFirst;
    int id=0; 
    
    stringstream tmp;
    tmp <<setfill('0')<<setw(2)<<id<<":"<<timerName;
    string timerNameID = tmp.str(); 
    HandleMap::const_iterator iter = handleMap_.find(timerNameID);
    if (iter != handleMap_.end())
      return  iter->second;
    else
      {
	int nThreads = omp_get_max_threads(); 
	for (int id =0;id<nThreads;id++)
	  {
	    stringstream tmp;
            tmp  <<setfill('0')<<setw(2)<<id<<":"<<timerName;
	    string timerNameID = tmp.str(); 
	    handle = timers_.size();
	    if (id==0) handleFirst=handle; 
	    handleMap_[timerNameID] = handle;
	    timers_.push_back(TimerStruct());
	  }
      }
    return handleFirst;
}

void profileSetRefTimer(const string& timerName)
{
   refTimer_ = timerName;
}

void profileSetPrintOrder(const string& timerName)
{
   printOrder_.push_back(timerName);
}



vector<string> generateOutputOrder(const vector<string>& words)
{
   vector<string> order;
   set<string> printed;
   
   int nThreads = omp_get_max_threads(); 
   for (int id =0;id<nThreads;id++)
   {
   stringstream tmp;
   tmp  <<setfill('0')<<setw(2)<<id<<":";
   string prefix = tmp.str(); 
   for (unsigned ii=0; ii<printOrder_.size(); ++ii)
   {
      string  name = prefix+printOrder_[ii]; 
	//printf("name=%s\n",name.c_str()); 
      if ( printOrder_[ii].empty() || 
         find(words.begin(), words.end(), name) != words.end() )
      {
         order.push_back(name);
         printed.insert(name);
      }
   }
   for (unsigned ii=0; ii<words.size(); ++ii)
   {
      if (printed.count(words[ii]) == 0 && prefix == words[ii].substr(0,3))
         order.push_back(words[ii]);
   }
   }
   //for (int ii=0;ii<order.size();ii++)  printf("order %d %s\n",ii,order[ii].c_str()); 
   return order;
}
vector<string> generateOutputOrder()
{
   vector<string> words;
   for (HandleMap::const_iterator iter=handleMap_.begin(); iter!=handleMap_.end(); ++iter)
   {
         words.push_back(iter->first);
   }
   //for (int ii=0;ii<words.size();ii++)  printf("%d %s\n",ii,words[ii].c_str()); 
   
   return generateOutputOrder(words);
}

void profileDumpTimes(ostream& out)
{
   double tick = getTick(); 
   string::size_type maxLen = 0;
   for (HandleMap::iterator iter=handleMap_.begin();
        iter!=handleMap_.end(); ++iter)
       maxLen = max(maxLen, iter->first.size());

   vector<string> outputOrder = generateOutputOrder();

   ios::fmtflags oldFlags = out.setf(ios::fixed, ios::floatfield);

   out << setw(maxLen+3) << " "
       << setw(10)     << "#Calls" << "  |"
       << setw(10)      << "Current" 
       << setw(10)      << "Average" 
       << setw(10)      << "Total"
       << "   "
       << setw(8) <<setprecision(2)     << "%Loop" << endl;
   out << "--------------------------------------------------------------------------------"
       << endl;

   double refTime = timers_[handleMap_[refTimer_]].total[CYCLES]*tick;

   for (unsigned iName=0; iName<outputOrder.size(); ++iName)
   {
      const string& name = outputOrder[iName];
      unsigned ii = handleMap_[name];
      if (timers_[ii].total[CYCLES]*tick < 0.0001) continue; 
      if (name.empty())
      {
         out << "--------------------------------------------------------------------------------" << endl;
         continue;
      }
      
      out << setw(maxLen) << left << name << " : "
          << setw(10) << right << timers_[ii].total[NCALLS] << "  |"
          << setw(10) << timers_[ii].total[CYCLES]*tick/timers_[ii].total[NCALLS] 
          << setw(10) << timers_[ii].total[CYCLES]*tick 
          << "   "
          << setw(8) << setprecision(2) << 100.0*((timers_[ii].total[CYCLES]*tick)/refTime)
          << endl;
   }
   out.setf(oldFlags);
   
}

void profileDumpAll(const string& dirname)
{
   MPI_Comm comm = MPI_COMM_WORLD;
   int myRank;
   MPI_Comm_rank(comm, &myRank);
   if (myRank == 0)
      DirTestCreate(dirname.c_str());
   string filename = dirname + "/profile";
   PFILE* file = Popen(filename.c_str(), "w", comm);

   Pprintf(file, "------------------------------\n");
   Pprintf(file, "Performance for rank %u\n", myRank);
   Pprintf(file, "------------------------------\n\n");
   stringstream buf;
   profileDumpTimes(buf);
   Pwrite(buf.str().c_str(),strlen(buf.str().c_str()),1,file);
   Pprintf(file, "\n\n");
   //printf("end profileDumpTimes\n"); 
   Pclose(file);
}

void profileDumpStats(ostream& out)
{
   MPI_Comm comm = MPI_COMM_WORLD;
   int myRank;
   MPI_Comm_rank(comm, &myRank);
   vector<string> timerName;
   for (HandleMap::const_iterator iter=handleMap_.begin();
        iter!=handleMap_.end(); ++iter)
      timerName.push_back(iter->first);
   unionOfStrings(timerName, comm, getUniqueTag(comm));

   vector<string> outputOrder = generateOutputOrder(timerName);
   unsigned nTimers = outputOrder.size();
   double tick = getTick(); 
   
   vector<double> aveTime(nTimers);
   vector<int>    nActive(nTimers);
   {
      unsigned bufSize = (1+2)*nTimers;
      uint64_t sendBuf[bufSize];
      
      for (unsigned ii=0; ii<nTimers; ++ii)
      {
         HandleMap::const_iterator here = handleMap_.find(outputOrder[ii]);
         if (here == handleMap_.end())
         {
            sendBuf[3*ii+0] = 0; // inactive timer
            sendBuf[3*ii+1] = 0;
            sendBuf[3*ii+2] = 0;
         }
         else
         {
	    uint64_t total=timers_[here->second].total[CYCLES];
            uint64_t lower=total&0xffffffff; 
            uint64_t upper=(total>>32)&0xffffffff; 
            sendBuf[3*ii+0] = 1; // active timer
            sendBuf[3*ii+1] = lower;
            sendBuf[3*ii+2] = upper;
         }
      }
      uint64_t recvBuf[bufSize];
      MPI_Allreduce(sendBuf, recvBuf, bufSize, MPI_LONG_LONG, MPI_SUM, comm);

      for (unsigned ii=0; ii<nTimers; ++ii)
      {
         nActive[ii] = (int) recvBuf[3*ii+0];
         uint64_t lower = recvBuf[3*ii+1]; 
         uint64_t upper = recvBuf[3*ii+2]; 
	 double time =  (upper*(tick*pow(2.0,32))) + (lower*tick);
         aveTime[ii] = (nActive[ii] > 0 ? (time/nActive[ii]) : 0.0);
      }
   }
   struct DoubleInt
   {
      double val;
      int rank;
   };
   
   DoubleInt minLoc[nTimers];
   DoubleInt maxLoc[nTimers];
   double sigma[nTimers];
   {
      DoubleInt minLocSendBuf[nTimers];
      DoubleInt maxLocSendBuf[nTimers];
      double sigmaSendBuf[nTimers];
      
      for (unsigned ii=0; ii<nTimers; ++ii)
      {
         minLocSendBuf[ii].rank = myRank;
         maxLocSendBuf[ii].rank = myRank;
         
         HandleMap::const_iterator here = handleMap_.find(outputOrder[ii]);
         if (here == handleMap_.end())
         {
            minLocSendBuf[ii].val = aveTime[ii];
            maxLocSendBuf[ii].val = aveTime[ii];
            sigmaSendBuf[ii] = 0;
         }
         else
         {
            //ewd: should these be set to zero for case when nActive[ii] == 0?
            minLocSendBuf[ii].val = timers_[here->second].total[CYCLES]*tick;
            maxLocSendBuf[ii].val = timers_[here->second].total[CYCLES]*tick;
            sigmaSendBuf[ii] = (timers_[here->second].total[CYCLES]*tick - aveTime[ii]);
            if (nActive[ii] > 0)
               sigmaSendBuf[ii] *= sigmaSendBuf[ii]/nActive[ii];
            else
               sigmaSendBuf[ii] = 0.0;
               
         }
      }
      MPI_Reduce(minLocSendBuf, minLoc, nTimers, MPI_DOUBLE_INT, MPI_MINLOC, 0, comm);
      MPI_Reduce(maxLocSendBuf, maxLoc, nTimers, MPI_DOUBLE_INT, MPI_MAXLOC, 0, comm);
      MPI_Reduce(sigmaSendBuf, sigma, nTimers, MPI_DOUBLE, MPI_SUM, 0, comm);
   }
   
   if (myRank != 0) return;

   string::size_type maxLen = 0;
   for (unsigned ii=0; ii<nTimers; ++ii)
      maxLen = max(maxLen, outputOrder[ii].size());
   ios::fmtflags oldFlags = out.setf(ios::fixed, ios::floatfield);
   out << setw(maxLen+3) << " "
       << setw(6) << "nTasks" << " |"
       << setw(10) << "minTime" << " ( rank )"
       << setw(10) << "maxTime" << " ( rank )"
       << setw(10) << "sigma"
       << setw(10) << "aveTime"
       << endl;
   out << "--------------------------------------------------------------------------------" << endl;
   for (unsigned ii=0; ii<nTimers; ++ii)
   {
      if (aveTime[ii] < 0.0001) continue;
      if (outputOrder[ii].empty())
      {
         out << "--------------------------------------------------------------------------------" << endl;
         continue;
      }
      out << setw(maxLen) << left << outputOrder[ii] << " : " << right
          << setw(6) << nActive[ii] << " |"
          << setprecision(3)
          << setw(10) << minLoc[ii].val
          << " (" << setw(6) << minLoc[ii].rank << ")"
          << setw(10) << maxLoc[ii].val
          << " (" << setw(6) << maxLoc[ii].rank << ")"
          << setw(10) << sqrt(sigma[ii])
          << setw(10) << aveTime[ii]
          << endl;
   }

   out.setf(oldFlags);
}
