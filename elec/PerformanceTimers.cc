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
#include <cstring>
#include <stdint.h>
#include "pio.h"
#include "mpiUtils.h"
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
   TimerHandle simulationLoopTimer;
   TimerHandle sensorTimer;
   TimerHandle sensorEvalTimer;
   TimerHandle sensorCompColorCenterTimer;
   TimerHandle sensorSetupLSTimer;
   TimerHandle sensorComputeLSTimer;
   TimerHandle sensorPrintTimer;
   TimerHandle haloTimer;
   TimerHandle haloLaunchTimer;
   TimerHandle haloWaitTimer;
   TimerHandle haloMove2BufTimer;
   TimerHandle diffusionCalcTimer;
   TimerHandle stimulusTimer;
   TimerHandle reactionTimer;
   TimerHandle reactionMiscTimer;
   TimerHandle nonGateTimer;
   TimerHandle nonGateRLTimer;
   TimerHandle GateNonGateTimer;
   TimerHandle gateTimer;
   TimerHandle gateRLTimer;
   TimerHandle diffusionLoopTimer;
   TimerHandle integratorTimer;
   TimerHandle reactionLoopTimer;
   TimerHandle reactionWaitTimer;
   TimerHandle diffusionWaitTimer;
   TimerHandle diffusionStallTimer;
   TimerHandle diffusionImbalanceTimer;
   TimerHandle imbalanceTimer;
   TimerHandle dummyTimer;
   TimerHandle diffusionL2BarrierHalo1Timer;
   TimerHandle diffusionL2BarrierHalo2Timer;
   TimerHandle diffusiondVmRCopyTimer;
   TimerHandle reactionL2ArriveTimer;
   TimerHandle reactionL2ResetTimer;
   TimerHandle FGR_ArrayLocal2MatrixTimer;
   TimerHandle FGR_ArrayRemote2MatrixTimer;
   TimerHandle FGR_BarrierTimer;
   TimerHandle FGR_Barrier2Timer;
   TimerHandle FGR_AlignCopyTimer;
   TimerHandle FGR_StencilTimer;
   TimerHandle FGR_Matrix2ArrayTimer;
   TimerHandle FGR_2D_StencilTimer;
   TimerHandle FGR_Boundary2MatrixTimer; 
   TimerHandle initializeDVmDTimer;
   TimerHandle rangeCheckTimer;
   TimerHandle barrier1Timer;
   TimerHandle barrier2Timer;
   TimerHandle barrier3Timer;
   TimerHandle barrier4Timer;
   TimerHandle barrier5Timer;
   TimerHandle printDataTimer;
   TimerHandle timingBarrierTimer;
   TimerHandle FGR_2D_StencilTimerTimer;
   TimerHandle stencilOverlapTimer;
   
   vector<TimerStruct> timers_;
   typedef map<string, TimerHandle> HandleMap;
   HandleMap handleMap_;
   vector<string> printOrder_;
   string refTimer_;
   bool allCounters_;
}

using namespace PerformanceTimers;

void  profileInit()
{
   allCounters_ = true;
   loopIOTimer = profileGetHandle("LoopIO");
   simulationLoopTimer = profileGetHandle("SimulationLoop");;
   sensorTimer = profileGetHandle("Sensors");
   sensorEvalTimer = profileGetHandle("SensorsEval");
   sensorCompColorCenterTimer = profileGetHandle("SensorCompColorCenter");
   sensorSetupLSTimer = profileGetHandle("SensorSetupLS");
   sensorComputeLSTimer = profileGetHandle("SensorComputeLS");
   sensorPrintTimer = profileGetHandle("SensorsPrint");
   haloTimer = profileGetHandle("HaloExchange");
   haloLaunchTimer = profileGetHandle("HaloExchangeLaunch");
   haloWaitTimer = profileGetHandle("HaloExchangeWait");
   haloMove2BufTimer = profileGetHandle("HaloExchMove2Buf");
   diffusionCalcTimer= profileGetHandle("DiffusionCalc");
   stimulusTimer = profileGetHandle("Stimulus");
   reactionTimer= profileGetHandle("Reaction");
   reactionMiscTimer= profileGetHandle("ReactionMisc");
   rangeCheckTimer= profileGetHandle("rangeCheck");
   nonGateTimer= profileGetHandle("Reaction_nonGate");
   nonGateRLTimer= profileGetHandle("Reaction_nonGateRL");
   gateTimer = profileGetHandle("Reaction_Gate");
   gateRLTimer = profileGetHandle("Reaction_GateRL");
   GateNonGateTimer = profileGetHandle("GateNonGateBarrier");
   diffusionLoopTimer= profileGetHandle("DiffusionLoop");
   integratorTimer = profileGetHandle("Integrator");
   haloTimer = profileGetHandle("HaloExchange");
   reactionLoopTimer = profileGetHandle("ReactionLoop");
   reactionWaitTimer = profileGetHandle("ReactionWait");
   diffusionWaitTimer = profileGetHandle("DiffusionWait");
   diffusionStallTimer = profileGetHandle("DiffusionStall");
   diffusionImbalanceTimer = profileGetHandle("DiffusionImbalance");
   imbalanceTimer = profileGetHandle("Imbalance");
   dummyTimer = profileGetHandle("Dummy");
   diffusionL2BarrierHalo1Timer = profileGetHandle("DiffL2BarrHalo1");
   diffusionL2BarrierHalo2Timer = profileGetHandle("DiffL2BarrHalo2");
   diffusiondVmRCopyTimer = profileGetHandle("DiffdVmReactionCopy");
   reactionL2ArriveTimer = profileGetHandle("reactionL2Arrive");
   reactionL2ResetTimer = profileGetHandle("reactionL2Reset");
   FGR_ArrayLocal2MatrixTimer = profileGetHandle("FGR_ALocal2Matrix");
   FGR_ArrayRemote2MatrixTimer = profileGetHandle("FGR_ARemote2Matrix");
   FGR_BarrierTimer = profileGetHandle("FGR_Barrier");
   FGR_Barrier2Timer = profileGetHandle("FGR_Barrier2");
   FGR_AlignCopyTimer = profileGetHandle("FGR_AlignCopy");
   FGR_StencilTimer = profileGetHandle("FGR_Stencil");
   FGR_Matrix2ArrayTimer = profileGetHandle("FGR_Matrix2Array");
   initializeDVmDTimer = profileGetHandle("zerodVmDArray");
   rangeCheckTimer = profileGetHandle("RangeCheck");
   barrier1Timer = profileGetHandle("Barrier1");
   barrier2Timer = profileGetHandle("Barrier2");
   barrier3Timer = profileGetHandle("Barrier3");
   barrier4Timer = profileGetHandle("Barrier4");
   barrier5Timer = profileGetHandle("Barrier5");
   printDataTimer = profileGetHandle("PrintData");
   timingBarrierTimer = profileGetHandle("TimingBarrier");
   FGR_2D_StencilTimer = profileGetHandle("FGR_2D_Stencil");
   stencilOverlapTimer = profileGetHandle("stencilOverlap");
   machineSpecficInit(); 
}
void profileStart(const TimerHandle& handle)
{
   int tid=omp_get_thread_num() ;
   int id=handle+tid; 
   timers_[id].start[CYCLES] = getTime();
   if (allCounters_)
   {
      for (int i=2; i<nCounters_;i++) 
      {
         uint64_t counter;
         readCounter(counterHandle[tid], i-2, &counter);   
         timers_[id].start[i] = counter; 
      }
   }
}

void profileStop(const TimerHandle& handle)
{
   int tid=omp_get_thread_num() ;
   int id=handle+tid; 
   timers_[id].total[NCALLS] += 1;
   uint64_t delta = getTime() - timers_[id].start[CYCLES];
   timers_[id].total[CYCLES] += delta;
   if (allCounters_)
   {
      for (int i=2; i<nCounters_;i++) 
      {
         uint64_t counter;
         readCounter(counterHandle[tid], i-2, &counter);
         uint64_t delta  = counter - timers_[id].start[i];
         timers_[id].total[i] += delta;
      }
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

void profileStart_HW(const TimerHandle& handle)
{
   allCounters_=1; 
   profileStart(handle);
   allCounters_=0; 
}

void profileStop_HW(const TimerHandle& handle)
{
   allCounters_=1; 
   profileStop(handle);
   allCounters_=0; 
}
void profileStart_HW(const std::string& timerName)
{
   allCounters_=1; 
   profileStart(profileGetHandle(timerName));
   allCounters_=0; 
}
void profileStop_HW(const std::string& timerName)
{
   allCounters_=1; 
   profileStop(profileGetHandle(timerName));
   allCounters_=0; 
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

void profileSetVerbosity(const bool verbosity)
{
   allCounters_ = verbosity;
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
         if (printOrder_[ii].empty())
            order.push_back("");
         
         string  name = prefix+printOrder_[ii]; 
         //printf("name=%s\n",name.c_str()); 
         if (find(words.begin(), words.end(), name) != words.end() )
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
   
   return generateOutputOrder(words);
}

void profileDumpTimes(ostream& out)
{
   char line[1024]; 
   double tick = getTick(); 
   string::size_type maxLen = 0;
   for (HandleMap::iterator iter=handleMap_.begin();
        iter!=handleMap_.end(); ++iter)
      maxLen = max(maxLen, iter->first.size());

   vector<string> outputOrder = generateOutputOrder();

   ios::fmtflags oldFlags = out.setf(ios::fixed, ios::floatfield);
   
   int counter=CYCLES; 
   {

      out << setw(12) << counterNames_[counter] 
          << setw(maxLen-9) << " "
          << setw(10)     << "#Calls" << "  |"
          << setw(10)      << "Average" 
          << setw(10)      << "Total"
          << "   "
          << setw(8) <<setprecision(2)     << "%Loop" << endl;
      out << "--------------------------------------------------------------------------------"
          << endl;

      double refCount = timers_[handleMap_[refTimer_]].total[CYCLES];

      for (unsigned iName=0; iName<outputOrder.size(); ++iName)
      {
         const string& name = outputOrder[iName];
         unsigned ii = handleMap_[name];
         //if (timers_[ii].total[counter]*tick < 0.0001) continue; 
         if (name.empty())
         {
            out << "--------------------------------------------------------------------------------" << endl;
            continue;
         }
      
         double  ncalls = (double)timers_[ii].total[NCALLS]; 
         double  count = (double)timers_[ii].total[counter]; 
         out << setw(maxLen) << left << name << " : "
             << setw(10) << right << ncalls << "  |"
             << setw(10) << timers_[ii].total[counter]*tick/ncalls 
             << setw(10) << timers_[ii].total[counter]*tick 
             << "   "
             << setw(8) << setprecision(2) << 100.0*((timers_[ii].total[counter])/refCount)
             << endl;
      }
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
   Pclose(file);
}


void profileDumpStats(ostream& out)
{
   MPI_Comm comm = MPI_COMM_WORLD;
   int myRank;
   MPI_Comm_rank(comm, &myRank);
   unsigned nTasks = getSize(0); 
   vector<string> timerName;
   for (HandleMap::const_iterator iter=handleMap_.begin();
        iter!=handleMap_.end(); ++iter)
      timerName.push_back(iter->first);
   unionOfStrings(timerName, comm, getUniqueTag(comm));

   vector<string> outputOrder = generateOutputOrder(timerName);
   unsigned nTimers = outputOrder.size();
   double tick = getTick(); 
   ios::fmtflags oldFlags = out.setf(ios::fixed, ios::floatfield);
   
   int nThreads=omp_get_max_threads(); 
   int nMarks = nTimers/nThreads; 
   int nCore =  getNCores(); 
   double refCount = timers_[handleMap_[refTimer_]].total[CYCLES];

   double aveCount[nTimers]; 
   double nActive[nTimers]; 
   double perfCount[nTimers][nCounters_];
   // ordinarily I'd declare these variables inside the scope of the
   // loop.  However, this seems to trigger a compiler bug on BGP.  So,
   // for now, they stay here.
   unsigned bufSize = 2*nTimers;
   long double sendBuf[bufSize];
   struct DoubleInt
   {
      double val;
      int rank;
   };
   
   DoubleInt minLoc[nTimers];
   DoubleInt maxLoc[nTimers];
   double sigma[nTimers];
   // end of variable that should be in loop scope
   int nPrintCounters_ = (allCounters_ ? nCounters_ : 2);
   nPrintCounters_=nCounters_; 
   for (int counter=0;counter<nPrintCounters_;counter++)
   {
      
      for (unsigned ii=0; ii<nTimers; ++ii)
      {
         HandleMap::const_iterator here = handleMap_.find(outputOrder[ii]);
         long double *buf = sendBuf+2*ii; 
         if (here == handleMap_.end())
         {
            buf[0] = 0; // inactive timer
            buf[1] = 0;
         }
         else
         {
            buf[0] = 1; // active timer
            buf[1] =timers_[here->second].total[counter];
         }
      }
      long double recvBuf[bufSize];
      MPI_Allreduce(sendBuf, recvBuf, bufSize, MPI_LONG_DOUBLE, MPI_SUM, comm);
      for (unsigned ii=0; ii<nTimers; ++ii)
      {
         double countSum = recvBuf[2*ii+1]; 
         perfCount[ii][counter] = (countSum+1e-100); 
       
      }
      if (counter!=CYCLES) continue; 
      for (unsigned ii=0; ii<nTimers; ++ii)
      {
         double count = sendBuf[2*ii+1]; 
         nActive[ii]= recvBuf[2*ii]; 
         double countSum = recvBuf[2*ii+1]; 
         aveCount[ii] = (nActive[ii] > 0 ? (countSum/nActive[ii]) : 0.0);
       
      }

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
               minLocSendBuf[ii].val = aveCount[ii];
               maxLocSendBuf[ii].val = aveCount[ii];
               sigmaSendBuf[ii] = 0;
            }
            else
            {
               //ewd: should these be set to zero for case when nActive[ii] == 0?
               minLocSendBuf[ii].val = timers_[here->second].total[counter];
               maxLocSendBuf[ii].val = timers_[here->second].total[counter];
               sigmaSendBuf[ii] =     (timers_[here->second].total[counter] - aveCount[ii]);
               if (nActive[ii] > 0)
                  sigmaSendBuf[ii] *= sigmaSendBuf[ii]/nActive[ii];
               else
                  sigmaSendBuf[ii] = 0.0;
               
            }
         }
         MPI_Reduce(minLocSendBuf, minLoc, nTimers, MPI_DOUBLE_INT, MPI_MINLOC, 0, comm);
         MPI_Reduce(maxLocSendBuf, maxLoc, nTimers, MPI_DOUBLE_INT, MPI_MAXLOC, 0, comm);
         MPI_Reduce(sigmaSendBuf, sigma, nTimers, MPI_DOUBLE, MPI_SUM, 0, comm);
   
         if (myRank == 0) 
         {
   

            string::size_type maxLen = 0;
            for (unsigned ii=0; ii<nTimers; ++ii)
               maxLen = max(maxLen, outputOrder[ii].size());
            out << setw(10) << counterNames_[counter] 
                << setw(maxLen-7) << " "
                << setw(6) << "nTasks" << " |"
                << setw(10) << "minTime" << " ( rank )"
                << setw(10) << "maxTime" << " ( rank )"
                << setw(10) << "sigma"
                << setw(10) << "aveTime"
                << endl;
            out << "--------------------------------------------------------------------------------------------------" << endl;
            for (unsigned ii=0; ii<nTimers; ++ii)
            {
               //if (aveCount[ii]*tick < 0.0001) continue;
               if (outputOrder[ii].empty())
               {
                  out << "--------------------------------------------------------------------------------------------------" << endl;
                  continue;
               }
               out << setw(maxLen) << left << outputOrder[ii] << " : " << right
                   << setw(6) << nActive[ii] << " |"
                   << setprecision(3)
                   << setw(10) << minLoc[ii].val*tick
                   << " (" << setw(6) << minLoc[ii].rank << ")"
                   << setw(10) << maxLoc[ii].val*tick
                   << " (" << setw(6) << maxLoc[ii].rank << ")"
                   << setw(10) << sqrt(sigma[ii])*tick
                   << setw(10) << aveCount[ii]*tick
                   << endl;
            }
         }
      }
   }
   if (myRank ==0 && nPrintCounters_ > 2) 
   {
      FILE *file=fopen("perfCount.data","w"); 
      int nMarks = nTimers/nThreads;
      int paraIndex = -1; 
      for (unsigned ii=0; ii<nTimers; ++ii) { if ( outputOrder[ii] == "00:parallelDiffReac" ) {paraIndex = ii; break; }}
      // This code doesn't work.  Commenting it out to avoid seeing a
      // completely bogus flog number.  - dfr
//       if (paraIndex != -1) 
//       {
//          double flopCount=0.0; 
//          for (unsigned ii=paraIndex; ii<nTimers; ii+=nMarks) flopCount += perfCount[ii][7]; 
//          double flops = flopCount / (perfCount[paraIndex][CYCLES]*tick);
//          printf("Flop Rate = %e GFLOP\n", flops*1e-9); 
//       }
      for (unsigned ii=0; ii<nTimers; ++ii)
      {
         int ompID=strtol((outputOrder[ii].substr(0,2)).c_str(),NULL,10);
             
         int coreID = ompID%nCore; 
         int timerID = ii%nMarks; 
         double cycle=0.0; 
         for (unsigned jj=0; jj<nTimers; ++jj)
         {
            int jjcoreID = ompID%nCore; 
            int jjtimerID = ii%nMarks; 
            if (jjcoreID == coreID && jjtimerID == timerID) cycle+=perfCount[jj][CYCLES]; 
         }
         double flop = 0.0 ; 
         if (cycle > 1e-9) flop = perfCount[ii][7]/(cycle*tick);
         
         if (perfCount[ii][NCALLS]>1e-9) 
         {
            fprintf(file,"%4d %-32s %2d %2d %3d %12.0f",ii,(outputOrder[ii].c_str())+3,ompID,coreID,timerID,perfCount[ii][0]/nTasks); 
            for (int jj=1;jj< nCounters_;jj++) fprintf(file," %7.3f", perfCount[ii][jj]*1e-9/nTasks); 
            fprintf(file," %13.6f", flop*1e-9); 
            fprintf(file,"\n"); 
            fflush(file); 
         }
      }
   }
   out.setf(oldFlags);
}
