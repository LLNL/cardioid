#include "PerformanceTimers.hh"
#include <vector>
#include <map>
#include <set>
#include <sys/time.h>
#include <iostream>
#include <iomanip>

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



vector<string> generatePrintOrder()
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

void profileDumpTimes(ostream& out)
{
   string::size_type maxLen = 0;
   for (HandleMap::iterator iter=handleMap_.begin();
        iter!=handleMap_.end(); ++iter)
       maxLen = max(maxLen, iter->first.size());

   vector<string> printOrder = generatePrintOrder();

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

   for (unsigned iName=0; iName<printOrder.size(); ++iName)
   {
      const string& name = printOrder[iName];
      if (name.empty())
         out << "--------------------------------------------------------------------------------" << endl;
      
      unsigned ii = handleMap_[name];
   out << setw(maxLen) << left << name << " : "
       << setw(10) << right << int(timers_[ii].nCalls) << "  |"
       << setw(10)  <<setprecision(3) << timers_[ii].lastTime 
       << setw(10)  << timers_[ii].totalTime/timers_[ii].nCalls 
       << setw(10)  << timers_[ii].totalTime
       << "   "
       << setw(8) << setprecision(2) << 100.0*(timers_[ii].totalTime/refTime)
       << endl;
   }
    out.setf(oldFlags);
    
}
