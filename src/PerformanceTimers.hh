#ifndef PERFORMANCE_TIMERS_HH
#define PERFORMANCE_TIMERS_HH

#include <string>
#include <iosfwd>

typedef unsigned TimerHandle;

namespace PerformanceTimers
{
   extern TimerHandle loopIOTimer;
   extern TimerHandle simulationLoopTimer;
   extern TimerHandle sensorTimer;
   extern TimerHandle sensorEvalTimer;
   extern TimerHandle sensorPrintTimer;
   extern TimerHandle sensorCompColorCenterTimer;
   extern TimerHandle sensorSetupLSTimer;
   extern TimerHandle sensorComputeLSTimer;
   extern TimerHandle haloTimer;
   extern TimerHandle haloLaunchTimer;
   extern TimerHandle haloWaitTimer;
   extern TimerHandle haloMove2BufTimer;
   extern TimerHandle diffusionCalcTimer;
   extern TimerHandle stimulusTimer;
   extern TimerHandle reactionTimer;
   extern TimerHandle reactionMiscTimer;
   extern TimerHandle nonGateTimer;
   extern TimerHandle nonGateRLTimer;
   extern TimerHandle GateNonGateTimer;
   extern TimerHandle gateTimer;
   extern TimerHandle gateRLTimer;
   extern TimerHandle diffusionLoopTimer;
   extern TimerHandle integratorTimer;
   extern TimerHandle reactionLoopTimer;
   extern TimerHandle reactionWaitTimer;
   extern TimerHandle diffusionWaitTimer;
   extern TimerHandle diffusionStallTimer;
   extern TimerHandle diffusionImbalanceTimer;
   extern TimerHandle imbalanceTimer;
   extern TimerHandle dummyTimer;
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
   extern TimerHandle FGR_2D_StencilTimer;
   extern TimerHandle FGR_Boundary2MatrixTimer; 
   extern TimerHandle initializeDVmDTimer;
   extern TimerHandle rangeCheckTimer;
   extern TimerHandle barrier1Timer;
   extern TimerHandle barrier2Timer;
   extern TimerHandle barrier3Timer;
   extern TimerHandle barrier4Timer;
   extern TimerHandle barrier5Timer;
   extern TimerHandle printDataTimer;
   extern TimerHandle timingBarrierTimer;
   extern TimerHandle stencilOverlapTimer;
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

void profileStart_HW(const TimerHandle& handle);
void profileStop_HW(const TimerHandle& handle);
void profileStart_HW(const std::string& timerName);
void profileStop_HW(const std::string& timerName);



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
