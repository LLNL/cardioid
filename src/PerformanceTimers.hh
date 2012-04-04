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
   extern TimerHandle diffusionCalcTimer;
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
   extern TimerHandle haloTimerExecute;
   extern TimerHandle haloTimerComplete;
   extern TimerHandle diffusionL2BarrierHalo1Timer;
   extern TimerHandle diffusionL2BarrierHalo2Timer;
   extern TimerHandle diffusiondVmRCopyTimer;
   extern TimerHandle reactionL2ArriveTimer;
   extern TimerHandle reactionL2ResetTimer;
   extern TimerHandle FGR_ArrayLocal2MatrixTimer;
   extern TimerHandle FGR_ArrayRemote2MatrixTimer;
   extern TimerHandle FGR_BarrierTimer;
   extern TimerHandle FGR_Barrier2Timer;
   extern TimerHandle FGR_AlignCopyTimer;
   extern TimerHandle FGR_StencilTimer;
   extern TimerHandle FGR_Matrix2ArrayTimer;
   extern TimerHandle haloMove2BufTimer;
};

/** Use the startTimer and stopTimer macros for timers that are inside
 * the main loop.  These can be compiled away by defining NTIMING.  For
 * timers outside the main loop consider calling profileStart and
 * profileStop instead */
#ifndef NTIMING
#define startTimer(handle)    \
   do                         \
   {                          \
      profileStart(handle);   \
   } while(0)
#define stopTimer(handle)     \
   do                         \
   {                          \
      profileStop(handle);    \
   } while(0)
#else
#define startTimer(handle)
#define stopTimer(handle)
#endif

/** Use these function calls only for timers that should *never* be
 *  turned off.  Typically this means they are outside the main
 *  simulation loop.  If the timer is inside the main loop use
 *  startTimer and stopTimer instead */
void profileStart(const TimerHandle& handle);
void profileStop(const TimerHandle& handle);
void profileStart(const std::string& timerName);
void profileStop(const std::string& timerName);




TimerHandle profileGetHandle(std::string timerName);

void profileSetVerbosity(const bool verbosity);
void profileSetRefTimer(const std::string& timerName);
void profileSetPrintOrder(const std::string& timerName);
void profileDumpTimes(std::ostream&);
void profileDumpAll(const std::string& dirname);
void profileDumpStats(std::ostream& out);
void profileInit();
void profilePrintName();
#endif
