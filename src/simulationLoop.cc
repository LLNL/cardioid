#include "simulationLoop.hh"

#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdio>
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
#include "threading.hh"

using namespace std;


static L2_Barrier_t reactionBarrier;
static L2_Barrier_t diffusionBarrier;
   
void simulationProlog(Simulate& sim)
{
   // initialize membrane voltage with a default value from the reaction model. 
   sim.VmArray_.resize(sim.anatomy_.size());
   sim.reaction_->initializeMembraneVoltage(sim.VmArray_);
  
   for (unsigned ii=sim.anatomy_.nLocal(); ii<sim.anatomy_.size(); ++ii) sim.VmArray_[ii] = 0;

   // Load state file, assign corresponding values to membrane voltage and cell model
   if (!sim.stateFilename_.empty()) readCheckpoint(sim, MPI_COMM_WORLD);
}
void loopIO(const Simulate& sim, const vector<double>& dVmR, vector<double>& dVmD,  const vector<double>& dVmE)
{
   int diffusionID = groupThreadID(sim.diffusionGroup_); 
   if (diffusionID != 0) return; 
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
      
   static TimerHandle loopIOHandle = profileGetHandle("LoopIO");
   static int firstCall=1; 
   int loop = sim.loop_; 

   // SENSORS
   static TimerHandle sensorHandle = profileGetHandle("Sensors");
   profileStart(sensorHandle);
   for (unsigned ii=0; ii<sim.sensor_.size(); ++ii)
   {
       if (loop % sim.sensor_[ii]->evalRate() == 0) sim.sensor_[ii]->eval(sim.time_, loop, sim.VmArray_, dVmR, dVmD, dVmE);
   }
   profileStop(sensorHandle);

   profileStart(loopIOHandle);
      
   for (unsigned ii=0; ii<sim.sensor_.size(); ++ii)
   {
      if (loop % sim.sensor_[ii]->printRate() == 0) sim.sensor_[ii]->print(sim.time_, loop, sim.VmArray_, dVmR, dVmD, dVmE);
   }
      
   if ( (loop % sim.printRate_ == 0) && myRank == 0)
   {
      if (firstCall) printf("    Loop     Time         Vm(t)        dVm_r(t-h)      dVm_d(t-h)       dVm_e(t-h)\n");
      printf("%8d %8.3f %15.8f %15.8f %15.8f %15.8f\n",loop,sim.time_,sim.VmArray_[0],dVmR[0],dVmD[0],dVmE[0]); 
   }
   if (!firstCall) 
   { 
      if (loop % sim.snapshotRate_ == 0)
      {
         stringstream name;
         name << "snapshot."<<setfill('0')<<setw(12)<<loop;
         string fullname = name.str();
         if (myRank == 0) DirTestCreate(fullname.c_str());
         fullname += "/anatomy";
         writeCells(sim, fullname.c_str());
      }
         
   }
   firstCall=0; 
   profileStop(loopIOHandle);
}


void simulationLoop(Simulate& sim)
{

   vector<double> dVmDiffusion(sim.anatomy_.nLocal(), 0.0);
   vector<double> dVmReaction(sim.anatomy_.nLocal(), 0.0);
   vector<double> dVmExternal(sim.anatomy_.nLocal(), 0.0);
   vector<double> iStim(sim.anatomy_.nLocal(), 0.0);
   simulationProlog(sim);
#ifdef SPI   
   spi_HaloExchange<double> voltageExchange(sim.sendMap_, (sim.commTable_));
#else
   mpi_HaloExchange<double> voltageExchange(sim.sendMap_, (sim.commTable_));
#endif
   while ( sim.loop_<sim.maxLoop_ )
   {
      int nLocal = sim.anatomy_.nLocal();
    
      static TimerHandle haloHandle = profileGetHandle("Halo Exchange");
      profileStart(haloHandle);
      voltageExchange.execute(sim.VmArray_, nLocal);
      voltageExchange.complete();
      profileStop(haloHandle);

      for (unsigned ii=0; ii<nLocal; ++ii) dVmExternal[ii] = 0;
    

      // DIFFUSION
      static TimerHandle diffusionHandle = profileGetHandle("Diffusion");
      profileStart(diffusionHandle);
      sim.diffusion_->calc(sim.VmArray_, dVmDiffusion, voltageExchange.get_recv_buf_(), nLocal);
      profileStop(diffusionHandle);

      // code to limit or set iStimArray goes here.
      static TimerHandle stimulusHandle = profileGetHandle("Stimulus");
      profileStart(stimulusHandle);
      for (unsigned ii=0; ii<sim.stimulus_.size(); ++ii) sim.stimulus_[ii]->stim(sim.time_, dVmDiffusion, dVmExternal);
      for (unsigned ii=0; ii<nLocal; ++ii) iStim[ii] = -(dVmDiffusion[ii] + dVmExternal[ii]);
      profileStop(stimulusHandle);

      
      // REACTION
      static TimerHandle reactionHandle = profileGetHandle("Reaction");
      profileStart(reactionHandle);
      sim.reaction_->calc(sim.dt_, sim.VmArray_, iStim, dVmReaction);
      profileStop(reactionHandle);

      static TimerHandle integratorHandle = profileGetHandle("Integrator");
      profileStart(integratorHandle);
      for (unsigned ii=0; ii<nLocal; ++ii)
      {
         double dVm = dVmReaction[ii] + dVmDiffusion[ii] + dVmExternal[ii] ;
         sim.VmArray_[ii] += sim.dt_*dVm;
      }
      sim.time_ += sim.dt_;
      ++sim.loop_;
      profileStop(integratorHandle);
      if (sim.checkpointRate_ > 0 && sim.loop_ % sim.checkpointRate_ == 0) writeCheckpoint(sim, MPI_COMM_WORLD);
      loopIO(sim,dVmReaction,dVmDiffusion,dVmExternal);

   }
}
void nullDiffusionLoop(Simulate& sim, vector<double>& dVmReaction, L2_BarrierHandle_t& reactionHandle, L2_BarrierHandle_t& diffusionHandle)
{
   int ompTid = omp_get_thread_num();
   //fprintf(tfile[ompTid],"nullDiffusionLoop: ompTid=%d\n",ompTid); fflush(stdout); 
   while ( sim.loop_<=sim.maxLoop_ )
   {
      L2_BarrierWithSync_WaitAndReset(&reactionBarrier, &reactionHandle, sim.reactionGroup_->nThreads);
      L2_BarrierWithSync_Arrive(&diffusionBarrier, &diffusionHandle, sim.diffusionGroup_->nThreads);
      L2_BarrierWithSync_Reset(&diffusionBarrier, &diffusionHandle, sim.diffusionGroup_->nThreads);
   }
}
void diffusionLoop(Simulate& sim, vector<double>& dVmReaction, L2_BarrierHandle_t& reactionHandle, L2_BarrierHandle_t& diffusionHandle)
{
    assert(groupThreadID(sim.diffusionGroup_)==0); //DiffusionLoop current only works on one thread. 
#ifdef SPI   
   spi_HaloExchange<double> voltageExchange(sim.sendMap_, (sim.commTable_));
#else
   mpi_HaloExchange<double> voltageExchange(sim.sendMap_, (sim.commTable_));
#endif
   int nLocal = sim.anatomy_.nLocal();
   vector<double> dVmExternal(nLocal, 0.0);
   vector<double> dVmDiffusion(nLocal, 0.0);
   vector<double> dVmReactionCpy(nLocal, 0.0);
   while ( sim.loop_ < sim.maxLoop_ )
   {
      int nLocal = sim.anatomy_.nLocal();
    
      static TimerHandle haloHandle = profileGetHandle("Halo Exchange");
      profileStart(haloHandle);
      voltageExchange.execute(sim.VmArray_, nLocal);
      voltageExchange.complete();
      profileStop(haloHandle);

      for (unsigned ii=0; ii<nLocal; ++ii) dVmDiffusion[ii] = 0;
      for (unsigned ii=0; ii<nLocal; ++ii) dVmExternal[ii] = 0;
    

      // DIFFUSION
      static TimerHandle diffusionTimerHandle = profileGetHandle("Diffusion");
      profileStart(diffusionTimerHandle);
      sim.diffusion_->calc(sim.VmArray_, dVmDiffusion, voltageExchange.get_recv_buf_(), nLocal);
      profileStop(diffusionTimerHandle);

      static TimerHandle stimulusHandle = profileGetHandle("Stimulus");
      profileStart(stimulusHandle);
      for (unsigned ii=0; ii<sim.stimulus_.size(); ++ii) sim.stimulus_[ii]->stim(sim.time_, dVmDiffusion, dVmExternal);
      profileStop(stimulusHandle);

      L2_BarrierWithSync_WaitAndReset(&reactionBarrier, &reactionHandle, sim.reactionGroup_->nThreads);

      dVmReactionCpy=dVmReaction; 
      static TimerHandle integratorHandle = profileGetHandle("Integrator");
      profileStart(integratorHandle);
      for (unsigned ii=0; ii<nLocal; ++ii)
      {
         double dVm = dVmReaction[ii] + dVmDiffusion[ii]+dVmExternal[ii];
         sim.VmArray_[ii] += sim.dt_*dVm;
      }
      sim.time_ += sim.dt_;
      ++sim.loop_;
      profileStop(integratorHandle);

      L2_BarrierWithSync_Arrive(&diffusionBarrier, &diffusionHandle, sim.diffusionGroup_->nThreads);
      L2_BarrierWithSync_Reset(&diffusionBarrier, &diffusionHandle, sim.diffusionGroup_->nThreads);
      if (sim.checkpointRate_ > 0 && sim.loop_ % sim.checkpointRate_ == 0) writeCheckpoint(sim, MPI_COMM_WORLD);
      loopIO(sim,dVmReactionCpy,dVmDiffusion,dVmExternal);
   }
}

void reactionLoop(Simulate& sim, vector<double>& dVmReaction,L2_BarrierHandle_t& reactionHandle, L2_BarrierHandle_t& diffusionHandle)
{
   while ( sim.loop_<sim.maxLoop_ )
   {
      int nLocal = sim.anatomy_.nLocal();
      
      //static TimerHandle reactionHandle = profileGetHandle("Reaction");
      //profileStart(reactionHandle);
      sim.reaction_->updateNonGate(sim.dt_, sim.VmArray_, dVmReaction);
      sim.reaction_->updateGate(sim.dt_, sim.VmArray_);
      // profileStop(reactionHandle);
      L2_BarrierWithSync_Arrive(&reactionBarrier, &reactionHandle, sim.reactionGroup_->nThreads);
      L2_BarrierWithSync_Reset(&reactionBarrier, &reactionHandle, sim.reactionGroup_->nThreads);
      L2_BarrierWithSync_WaitAndReset(&diffusionBarrier, &diffusionHandle, sim.diffusionGroup_->nThreads);
   }
}
void nullReactionLoop(Simulate& sim, vector<double>& dVmReaction,L2_BarrierHandle_t& reactionHandle, L2_BarrierHandle_t& diffusionHandle)
{
   while ( sim.loop_<=sim.maxLoop_ )
   {
      
      L2_BarrierWithSync_Arrive(&reactionBarrier, &reactionHandle, sim.reactionGroup_->nThreads);
      L2_BarrierWithSync_Reset(&reactionBarrier, &reactionHandle, sim.reactionGroup_->nThreads);
      L2_BarrierWithSync_WaitAndReset(&diffusionBarrier, &diffusionHandle, sim.diffusionGroup_->nThreads);
   }
}


void simulationLoopParallelDiffusionReaction(Simulate& sim)
{
   int nLocal = sim.anatomy_.nLocal();
   vector<double> dVmReaction(nLocal, 0.0);
   simulationProlog(sim);
#pragma omp parallel
   {
      int ompTid = omp_get_thread_num();
      
      
      L2_BarrierHandle_t reactionHandle;
      L2_BarrierHandle_t diffusionHandle;
      L2_BarrierWithSync_InitInThread(&reactionBarrier, &reactionHandle);
      L2_BarrierWithSync_InitInThread(&diffusionBarrier, &diffusionHandle);
      
#pragma omp barrier
      if ( sim.tinfo_.threadingMap_[ompTid] == sim.diffusionGroup_) 
      {
          int dID = groupThreadID(sim.diffusionGroup_); 
          if (dID ==0) diffusionLoop(sim, dVmReaction,  reactionHandle, diffusionHandle);
          else nullDiffusionLoop(sim, dVmReaction,  reactionHandle, diffusionHandle);
      }
      if ( sim.tinfo_.threadingMap_[ompTid] == sim.reactionGroup_) 
      {
          int rID = groupThreadID(sim.reactionGroup_); 
          reactionLoop(sim, dVmReaction, reactionHandle, diffusionHandle);
      } 
   }
}
