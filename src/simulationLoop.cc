#include "simulationLoop.hh"

#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <omp.h>

#include "Simulate.hh"
#include "Diffusion.hh"
#include "Reaction.hh"
#include "Stimulus.hh"
#include "Sensor.hh"
#include "HaloExchange.hh"
#include "GridRouter.hh"
#include "ioUtils.h"
#include "writeCells.hh"
#include "checkpointIO.hh"
#include "PerformanceTimers.hh"
#include "fastBarrier.hh"
#include "object_cc.hh"
#include "clooper.h"
#include "ThreadUtils.hh"


using namespace std;
using namespace PerformanceTimers;


void simulationProlog(Simulate& sim)
{
   // initialize membrane voltage with default value from the reaction model. 
   sim.VmArray_.resize(sim.anatomy_.size());
   sim.reaction_->initializeMembraneVoltage(sim.VmArray_);
  
   for (unsigned ii=sim.anatomy_.nLocal(); ii<sim.anatomy_.size(); ++ii)
      sim.VmArray_[ii] = 0;

   // Load state file, assign corresponding values to membrane voltage and cell model
   for (unsigned ii=0; ii<sim.stateFilename_.size(); ++ii)
      readCheckpoint(sim.stateFilename_[ii], sim, MPI_COMM_WORLD);

}


void loopIO(const Simulate& sim, const vector<double>& dVmR, vector<double>& dVmD)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
      
   static int firstCall=1; 
   int loop = sim.loop_; 

   // SENSORS
   startTimer(sensorTimer);
   for (unsigned ii=0; ii<sim.sensor_.size(); ++ii)
   {
       if (loop % sim.sensor_[ii]->evalRate() == 0) sim.sensor_[ii]->eval(sim.time_, loop, sim.VmArray_, dVmR, dVmD);
   }
   stopTimer(sensorTimer);

   startTimer(loopIOTimer);
   for (unsigned ii=0; ii<sim.sensor_.size(); ++ii)
   {
      if (loop % sim.sensor_[ii]->printRate() == 0) sim.sensor_[ii]->print(sim.time_, loop, sim.VmArray_, dVmR, dVmD);
   }
      
   if ( (loop % sim.printRate_ == 0) && myRank == 0)
   {
      if (firstCall) printf("    Loop     Time         Vm(t)        dVm_r(t-h)      dVm_d(t-h)\n");
      printf("%8d %8.3f %15.8f %15.8f %15.8f\n",loop,sim.time_,sim.VmArray_[0],dVmR[0],dVmD[0]); 
   }

   if (sim.checkpointRate_ > 0 && sim.loop_ % sim.checkpointRate_ == 0)
      writeCheckpoint(sim, MPI_COMM_WORLD);

   if (!firstCall) 
   { 
      if (loop % sim.snapshotRate_ == 0)
      {
         stringstream name;
         name << "snapshot."<<setfill('0')<<setw(12)<<loop;
         string fullname = name.str();
         if (myRank == 0) DirTestCreate(fullname.c_str());
         fullname += "/anatomy";
         if( sim.coarsedata_==0 )
            writeCells(sim, fullname.c_str());
         else{
            sim.coarsedata_->computeColorAverages(sim.VmArray_);
            sim.coarsedata_->writeAverages(fullname,sim.time_, sim.loop_);
         }
      }
         
   }
   firstCall=0; 
   stopTimer(loopIOTimer);
}


void simulationLoop(Simulate& sim)
{

   vector<double> dVmDiffusion(sim.anatomy_.nLocal(), 0.0);
   vector<double> dVmReaction(sim.anatomy_.nLocal(), 0.0);
   vector<double> iStim(sim.anatomy_.nLocal(), 0.0);
   simulationProlog(sim);
   HaloExchange<double> voltageExchange(sim.sendMap_, (sim.commTable_));

#if defined(SPI) && defined(TRACESPI)
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   cout << "Rank[" << myRank << "]: numOfNeighborToSend=" << sim.commTable_->_sendTask.size() << " numOfNeighborToRecv=" << sim.commTable_->_recvTask.size() << " numOfBytesToSend=" << sim.commTable_->_sendOffset[sim.commTable_->_sendTask.size()]*sizeof(double) << " numOfBytesToRecv=" << sim.commTable_->_recvOffset[sim.commTable_->_recvTask.size()]*sizeof(double) << endl;
#endif

   while ( sim.loop_<sim.maxLoop_ )
   {
      int nLocal = sim.anatomy_.nLocal();
    
      startTimer(imbalanceTimer);
      voltageExchange.barrier();
      stopTimer(imbalanceTimer);
      
      startTimer(haloTimer);
      voltageExchange.fillSendBuffer(sim.VmArray_);
      voltageExchange.startComm();
      voltageExchange.wait();
      stopTimer(haloTimer);

      // DIFFUSION
      startTimer(diffusionCalcTimer);
      sim.diffusion_->updateLocalVoltage(&(sim.VmArray_[0]));
      sim.diffusion_->updateRemoteVoltage(voltageExchange.getRecvBuf());
      sim.diffusion_->calc(dVmDiffusion);
      stopTimer(diffusionCalcTimer);

      startTimer(stimulusTimer);
      // add stimulus to dVmDiffusion
      for (unsigned ii=0; ii<sim.stimulus_.size(); ++ii)
         sim.stimulus_[ii]->stim(sim.time_, dVmDiffusion);
      for (unsigned ii=0; ii<nLocal; ++ii)
         iStim[ii] = -(dVmDiffusion[ii]);
      stopTimer(stimulusTimer);

      
      // REACTION
      startTimer(reactionTimer);
      sim.reaction_->calc(sim.dt_, sim.VmArray_, iStim, dVmReaction);
      stopTimer(reactionTimer);

      startTimer(integratorTimer);
      // no special BGQ integrator is this loop.  Mare bang for buck
      // from OMP threading.
      #pragma omp parallel for
      for (int ii=0; ii<nLocal; ++ii)
      {
         double dVm = dVmReaction[ii] + dVmDiffusion[ii];
         sim.VmArray_[ii] += sim.dt_*dVm;
      }

      sim.time_ += sim.dt_;
      ++sim.loop_;
      stopTimer(integratorTimer);
      loopIO(sim,dVmReaction,dVmDiffusion);
   }
}


/** One stop shopping for all of the data that we would rather create
 *  before the parallel section, either because we want the data to
 *  be shared across threads, or because it is just easier to create
 *  it in a single thread environment.  This data is all gathered
 *  into a struct so that it can be easily passed to the diffusion
 *  and reaction loops without giant argument lists.
 *
 *  You might argue that all of this data could be just as well
 *  stored in the Simulation class.  We choose not to that because
 *  - It would turn the Simulation class into quite a kitchen sink
 *  - It would make it even harder to construct the Simulation class
 *    properly
 *  - The data is only needed in the simulation loop.
 */
struct SimLoopData
{
   SimLoopData(const Simulate& sim)
   : voltageExchange(sim.sendMap_, (sim.commTable_))
   {
      haloBarrier = L2_BarrierWithSync_InitShared();
      reactionBarrier = L2_BarrierWithSync_InitShared();
      diffusionBarrier= L2_BarrierWithSync_InitShared();
      reactionWaitOnNonGateBarrier= L2_BarrierWithSync_InitShared();
      int nLocal = sim.anatomy_.nLocal();
      dVmDiffusion.resize(nLocal, 0.0);
      dVmReactionCpy.resize(nLocal, 0.0);
      dVmReaction.resize(nLocal, 0.0); 
      mkOffsets(integratorOffset, sim.anatomy_.nLocal(), sim.reactionThreads_);
   }

   ~SimLoopData()
   {
      free(reactionWaitOnNonGateBarrier);
      free(diffusionBarrier);
      free(reactionBarrier);
   }
   
   
   L2_Barrier_t* haloBarrier;
   L2_Barrier_t* reactionBarrier;
   L2_Barrier_t* diffusionBarrier;
   L2_Barrier_t* reactionWaitOnNonGateBarrier;
   
   vector<double> dVmReaction; 
   vector<double> dVmDiffusion;
   vector<int> integratorOffset;
   vector<double> dVmReactionCpy;
   HaloExchange<double> voltageExchange;
};

void diffusionLoop(Simulate& sim,
                   SimLoopData& loopData,
                   L2_BarrierHandle_t& reactionHandle,
                   L2_BarrierHandle_t& diffusionHandle)
{
   profileStart(diffusionLoopTimer);
   int tid = sim.diffusionThreads_.teamRank();
   L2_BarrierHandle_t haloBarrierHandle;
   L2_BarrierWithSync_InitInThread(loopData.haloBarrier, &haloBarrierHandle);

   while ( sim.loop_ < sim.maxLoop_ )
   {
      int nLocal = sim.anatomy_.nLocal();
      
      // Halo Exchange
      if (tid == 0)
      {
         startTimer(diffusionImbalanceTimer);
         loopData.voltageExchange.barrier();
         stopTimer(diffusionImbalanceTimer);

         startTimer(haloTimer);

         startTimer(haloTimerExecute);
         loopData.voltageExchange.fillSendBuffer(sim.VmArray_);
         loopData.voltageExchange.startComm();
         stopTimer(haloTimerExecute);

         startTimer(haloTimerComplete);
         loopData.voltageExchange.wait();
         stopTimer(haloTimerComplete);

         stopTimer(haloTimer);
      }
      
      // Need a barrier for the completion of the halo exchange.
      startTimer(diffusionL2BarrierHalo1Timer);
      L2_BarrierWithSync_Barrier(loopData.haloBarrier, &haloBarrierHandle,
                                 sim.diffusionThreads_.nThreads());
      stopTimer(diffusionL2BarrierHalo1Timer);
      
      startTimer(diffusionCalcTimer);
      // copy remote
      sim.diffusion_->updateRemoteVoltage(loopData.voltageExchange.getRecvBuf());
      // temporary barrier
      L2_BarrierWithSync_Barrier(loopData.haloBarrier, &haloBarrierHandle,
                                 sim.diffusionThreads_.nThreads());
      //stencil
      sim.diffusion_->calc(loopData.dVmDiffusion);
      stopTimer(diffusionCalcTimer);
      
      //barrier
      startTimer(diffusionL2BarrierHalo2Timer);
      L2_BarrierWithSync_Barrier(loopData.haloBarrier, &haloBarrierHandle,
                                 sim.diffusionThreads_.nThreads());
      stopTimer(diffusionL2BarrierHalo2Timer);

      // stimulus
      if (tid == 0)
      {
         startTimer(stimulusTimer);
         for (unsigned ii=0; ii<sim.stimulus_.size(); ++ii)
            sim.stimulus_[ii]->stim(sim.time_, loopData.dVmDiffusion);
         stopTimer(stimulusTimer);
      }

      // announce that diffusion derivatives are ready.
      L2_BarrierWithSync_Arrive(loopData.diffusionBarrier, &diffusionHandle,
                                sim.diffusionThreads_.nThreads());
      L2_BarrierWithSync_Reset(loopData.diffusionBarrier,
                               &diffusionHandle,
                               sim.diffusionThreads_.nThreads());

      
      // wait for reaction (integration) to finish
      startTimer(diffusionWaitTimer); 
      L2_BarrierWithSync_WaitAndReset(loopData.reactionBarrier,
                                      &reactionHandle,
                                      sim.reactionThreads_.nThreads());
      stopTimer(diffusionWaitTimer); 

      startTimer(dummyTimer);
      stopTimer(dummyTimer);

      if (tid == 0)
      {
         startTimer(diffusiondVmRCopyTimer);
         loopData.dVmReactionCpy = loopData.dVmReaction; 
         stopTimer(diffusiondVmRCopyTimer);

         loopIO(sim, loopData.dVmReactionCpy, loopData.dVmDiffusion);
      }
   }
   profileStop(diffusionLoopTimer);
}

void reactionLoop(Simulate& sim, SimLoopData& loopData, L2_BarrierHandle_t& reactionHandle, L2_BarrierHandle_t& diffusionHandle)
{
   profileStart(reactionLoopTimer);
   int tid = sim.reactionThreads_.teamRank();
   vector<double>& dVmReaction = loopData.dVmReaction;
   L2_BarrierHandle_t reactionWaitOnNonGateHandle;
   L2_BarrierWithSync_InitInThread(loopData.reactionWaitOnNonGateBarrier, &reactionWaitOnNonGateHandle);

   while ( sim.loop_<sim.maxLoop_ )
   {
      int nLocal = sim.anatomy_.nLocal();
      
      startTimer(reactionTimer);
      sim.reaction_->updateNonGate(sim.dt_, sim.VmArray_, dVmReaction);
      L2_BarrierWithSync_Barrier(loopData.reactionWaitOnNonGateBarrier, &reactionWaitOnNonGateHandle, sim.reactionThreads_.nThreads());
      sim.reaction_->updateGate(sim.dt_, sim.VmArray_);
      L2_BarrierWithSync_Barrier(loopData.reactionWaitOnNonGateBarrier, &reactionWaitOnNonGateHandle, sim.reactionThreads_.nThreads());
      stopTimer(reactionTimer);


      startTimer(dummyTimer);
      stopTimer(dummyTimer);

      startTimer(reactionWaitTimer);
      L2_BarrierWithSync_WaitAndReset(loopData.diffusionBarrier, &diffusionHandle, sim.diffusionThreads_.nThreads());
      stopTimer(reactionWaitTimer);

      startTimer(integratorTimer);
      integrateLoop(loopData.integratorOffset[tid],
                    loopData.integratorOffset[tid+1], sim.dt_,
                    &loopData.dVmReaction[0], &loopData.dVmDiffusion[0], &sim.VmArray_[0]);
      
      if (tid == 0)
      {
         sim.time_ += sim.dt_;
         ++sim.loop_;
      }

      stopTimer(integratorTimer);
      sim.diffusion_->updateLocalVoltage(&sim.VmArray_[0]);

      startTimer(reactionL2ArriveTimer);
      L2_BarrierWithSync_Arrive(loopData.reactionBarrier, &reactionHandle, sim.reactionThreads_.nThreads());
      stopTimer(reactionL2ArriveTimer);

      startTimer(reactionL2ResetTimer);
      L2_BarrierWithSync_Reset(loopData.reactionBarrier, &reactionHandle, sim.reactionThreads_.nThreads());
      stopTimer(reactionL2ResetTimer);
   }
   profileStop(reactionLoopTimer);
}

void nullReactionLoop(Simulate& sim, SimLoopData& loopData, L2_BarrierHandle_t& reactionHandle, L2_BarrierHandle_t& diffusionHandle)
{
   while ( sim.loop_<=sim.maxLoop_ )
   {
      L2_BarrierWithSync_Arrive(loopData.reactionBarrier, &reactionHandle, sim.reactionThreads_.nThreads());
      L2_BarrierWithSync_Reset(loopData.reactionBarrier, &reactionHandle, sim.reactionThreads_.nThreads());
      L2_BarrierWithSync_WaitAndReset(loopData.diffusionBarrier, &diffusionHandle, sim.diffusionThreads_.nThreads());
   }
}



void simulationLoopParallelDiffusionReaction(Simulate& sim)
{
   SimLoopData loopData(sim);

   simulationProlog(sim);

#if defined(SPI) && defined(TRACESPI)
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   cout << "Rank[" << myRank << "]: numOfNeighborToSend=" << sim.commTable_->_sendTask.size() << " numOfNeighborToRecv=" << sim.commTable_->_recvTask.size() << " numOfBytesToSend=" << sim.commTable_->_sendOffset[sim.commTable_->_sendTask.size()]*sizeof(double) << " numOfBytesToRecv=" << sim.commTable_->_recvOffset[sim.commTable_->_recvTask.size()]*sizeof(double) << endl;
#endif

   #pragma omp parallel
   {
      int ompTid = omp_get_thread_num();
      
      L2_BarrierHandle_t reactionHandle;
      L2_BarrierHandle_t diffusionHandle;
      L2_BarrierWithSync_InitInThread(loopData.reactionBarrier, &reactionHandle);
      L2_BarrierWithSync_InitInThread(loopData.diffusionBarrier, &diffusionHandle);
      
      // setup matrix voltages for first timestep.
      sim.diffusion_->updateLocalVoltage(&sim.VmArray_[0]);
      #pragma omp barrier
      profileStart(parallelDiffReacTimer);
      if ( sim.diffusionThreads_.teamRank() >= 0) 
      {
         diffusionLoop(sim, loopData, reactionHandle, diffusionHandle);
      }
      if ( sim.reactionThreads_.teamRank() >= 0) 
      {
          reactionLoop(sim, loopData, reactionHandle, diffusionHandle);
      } 
      profileStop(parallelDiffReacTimer);
   }
}
