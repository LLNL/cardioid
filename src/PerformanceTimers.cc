#include "PerformanceTimers.hh"
#include <vector>
#include <map>
#include <set>
#include <sys/time.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include "pio.h"
#include "ioUtils.h"
#include "tagServer.h"
#include "unionOfStrings.hh"

using namespace std;

namespace PerformanceTimers
{
   struct TimerStruct
   {
      double startTime;
      double totalTime;
      double lastTime;
      double nCalls;
   };
   
   vector<TimerStruct> timers_;
   typedef map<string, TimerHandle> HandleMap;
   HandleMap handleMap_;
   vector<string> printOrder_;
   string refTimer_;

   
   double getTime()
   {
      double t;
#ifdef BGL
      unsigned long long rts_get_timebase(void);
      double seconds_per_cycle = 1.4285714285714285714e-9;
      t = ((double)rts_get_timebase())*seconds_per_cycle;
#else
      struct timeval ptime;
      gettimeofday(&ptime, (struct timezone *)NULL);
      t = ptime.tv_sec + 1e-6*ptime.tv_usec;
#endif
      return t;
   }
}

using namespace PerformanceTimers;


       
void profileStart(const TimerHandle& handle)
{
   timers_[handle].startTime = getTime();
   ++timers_[handle].nCalls;
}

void profileStop(const TimerHandle& handle)
{
   timers_[handle].lastTime = getTime() - timers_[handle].startTime;
   timers_[handle].totalTime += timers_[handle].lastTime;
}

void profileStart(const std::string& timerName)
{
   profileStart(profileGetHandle(timerName));
}

void profileStop(const std::string& timerName)
{
   profileStop(profileGetHandle(timerName));
}

TimerHandle profileGetHandle(const string& timerName)
{
   HandleMap::const_iterator iter = handleMap_.find(timerName);
   if (iter != handleMap_.end())
      return iter->second;

   TimerHandle handle = timers_.size();
   handleMap_[timerName] = handle;
   timers_.push_back(TimerStruct());
   return handle;
}

void profileSetRefTimer(const string& timerName)
{
   refTimer_ = timerName;
}

void profileSetPrintOrder(const string& timerName)
{
   printOrder_.push_back(timerName);
}



vector<string> generateOutputOrder()
{
   vector<string> order;
   set<string> printed;
   
   for (unsigned ii=0; ii<printOrder_.size(); ++ii)
      if ( printOrder_[ii].empty() ||
           handleMap_.count(printOrder_[ii]) == 1 )
      {
         order.push_back(printOrder_[ii]);
         printed.insert(printOrder_[ii]);
      }
   
   for (HandleMap::const_iterator iter=handleMap_.begin();
        iter!=handleMap_.end(); ++iter)
   {
      if (printed.count(iter->first) == 0)
         order.push_back(iter->first);
   }
   
   return order;
}

vector<string> generateOutputOrder(const vector<string>& words)
{
   vector<string> order;
   set<string> printed;
   
   for (unsigned ii=0; ii<printOrder_.size(); ++ii)
      if ( printOrder_[ii].empty() ||
           find(words.begin(), words.end(), printOrder_[ii]) != words.end() )
      {
         order.push_back(printOrder_[ii]);
         printed.insert(printOrder_[ii]);
      }
   
   for (unsigned ii=0; ii<words.size(); ++ii)
   {
      if (printed.count(words[ii]) == 0)
         order.push_back(words[ii]);
   }
   
   return order;
}

void profileDumpTimes(ostream& out)
{
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

   double refTime = timers_[handleMap_[refTimer_]].totalTime;

   for (unsigned iName=0; iName<outputOrder.size(); ++iName)
   {
      const string& name = outputOrder[iName];
      if (name.empty())
      {
         out << "--------------------------------------------------------------------------------" << endl;
         continue;
      }
      
      unsigned ii = handleMap_[name];
      out << setw(maxLen) << left << name << " : "
          << setw(10) << right << int(timers_[ii].nCalls) << "  |"
          << setw(10) << setprecision(3) << timers_[ii].lastTime 
          << setw(10) << timers_[ii].totalTime/timers_[ii].nCalls 
          << setw(10) << timers_[ii].totalTime
          << "   "
          << setw(8) << setprecision(2) << 100.0*(timers_[ii].totalTime/refTime)
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
   Pprintf(file, "%s", buf.str().c_str());
   Pprintf(file, "\n\n");
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

   vector<double> aveTime(nTimers);
   vector<int>    nActive(nTimers);
   {
      unsigned bufSize = 2*nTimers;
      double sendBuf[bufSize];
      
      for (unsigned ii=0; ii<nTimers; ++ii)
      {
         HandleMap::const_iterator here = handleMap_.find(outputOrder[ii]);
         if (here == handleMap_.end())
         {
            sendBuf[2*ii+0] = 0; // inactive timer
            sendBuf[2*ii+1] = 0;
         }
         else
         {
            sendBuf[2*ii+0] = 1; // active timer
            sendBuf[2*ii+1] = timers_[here->second].totalTime;
         }
      }
      double recvBuf[bufSize];
      MPI_Allreduce(sendBuf, recvBuf, bufSize, MPI_DOUBLE, MPI_SUM, comm);
      for (unsigned ii=0; ii<nTimers; ++ii)
      {
         nActive[ii] = (int) recvBuf[2*ii+0];
         aveTime[ii] = recvBuf[2*ii+1]/nActive[ii];
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
            minLocSendBuf[ii].val = timers_[here->second].totalTime;
            maxLocSendBuf[ii].val = timers_[here->second].totalTime;
            sigmaSendBuf[ii] = (timers_[here->second].totalTime - aveTime[ii]);
            sigmaSendBuf[ii] *= sigmaSendBuf[ii]/nActive[ii];
         }
      }
      MPI_Reduce(minLocSendBuf, minLoc, nTimers, MPI_DOUBLE_INT, MPI_MINLOC, 0, comm);
      MPI_Reduce(maxLocSendBuf, maxLoc, nTimers, MPI_DOUBLE_INT, MPI_MAXLOC, 0, comm);
      MPI_Reduce(sigmaSendBuf, sigma, nTimers, MPI_DOUBLE, MPI_SUM, 0, comm);
   }
   
   if (myRank != 0)
      return;

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
