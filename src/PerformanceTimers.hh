#ifndef PERFORMANCE_TIMERS_HH
#define PERFORMANCE_TIMERS_HH

#include <string>
#include <iosfwd>

typedef unsigned TimerHandle;

namespace PerformanceTimers
{
   extern TimerHandle loopIOTimer;
   extern TimerHandle sensorTimer;
   extern TimerHandle haloTimer;
   extern TimerHandle diffusionTimer;
   extern TimerHandle stimulusTimer;
   extern TimerHandle reactionTimer;
   extern TimerHandle nonGateTimer;
   extern TimerHandle gateTimer;
   extern TimerHandle diffusionLoopTimer;
   extern TimerHandle integratorTimer;
   extern TimerHandle reactionLoopTimer;
   extern TimerHandle reactionWaitTimer;
   extern TimerHandle diffusionWaitTimer;
   extern TimerHandle diffusionStallTimer;
   extern TimerHandle diffusionImbalanceTimer;
   extern TimerHandle imbalanceTimer;
   extern TimerHandle dummyTimer;
   extern TimerHandle parallelDiffReacTimer;
};

void profileStart(const TimerHandle& handle);
void profileStop(const TimerHandle& handle);
void profileStart(const std::string& timerName);
void profileStop(const std::string& timerName);

TimerHandle profileGetHandle(std::string timerName);

void profileSetRefTimer(const std::string& timerName);
void profileSetPrintOrder(const std::string& timerName);
void profileDumpTimes(std::ostream&);
void profileDumpAll(const std::string& dirname);
void profileDumpStats(std::ostream& out);
void profileInit();
void profilePrintName();
#endif
