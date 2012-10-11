#include "simulationLoop.hh"

#include <vector>
#include <utility>
#include <algorithm>
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
#include "mpiUtils.h"

/*
  This follwing header, fastBarrier_nosync.hh, contains barrier code
  stolen from fastBarrierBGQ.hh, but with memory synchronization
  turned off. It is acceptable to use that when a barrier executes
  entirely within one core. It is necessary that a squad runs within
  one core.
  
  fastBarrier_nosync also sets the PER_SQUAD_BARRIER macro. By not
  including it, the original code with a barrier accross all reaction
  threads is used.
*/
//#include "fastBarrier_nosync.hh"

#include "object_cc.hh"
#include "clooper.h"
#include "ThreadUtils.hh"

using namespace std;
using namespace PerformanceTimers;

namespace
{
   void threadBarrier(TimerHandle tHandle, L2_Barrier_t* bPtr, L2_BarrierHandle_t* bHandlePtr, int nThreads)
   {
      startTimer(tHandle);
      L2_BarrierWithSync_Barrier(bPtr, bHandlePtr, nThreads);
      stopTimer(tHandle);
   }
}


void simulationProlog(Simulate& sim)
{
   // initialize membrane voltage with default value from the reaction
   // model.  Initialize with zero for remote cells
   sim.vdata_.setup(sim.anatomy_);
   sim.reaction_->initializeMembraneVoltage(sim.vdata_.VmArray_);
   for (unsigned ii=sim.anatomy_.nLocal(); ii<sim.anatomy_.size(); ++ii)
      sim.vdata_.VmArray_[ii] = 0;

   // Load state file, assign corresponding values to membrane voltage
   // and cell model.  This may overwrite the initialization we just
   // did.
   for (unsigned ii=0; ii<sim.stateFilename_.size(); ++ii)
      readCheckpoint(sim.stateFilename_[ii], sim, MPI_COMM_WORLD);
}


void loopIO(const Simulate& sim,int firstCall)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   const Anatomy& anatomy = sim.anatomy_;
      
   int loop = sim.loop_; 

   const VectorDouble32& VmArray(sim.vdata_.VmArray_);
   const VectorDouble32& dVmR(sim.vdata_.dVmReaction_);
   const VectorDouble32& dVmD(sim.vdata_.dVmDiffusion_);

   // SENSORS
#pragma omp critical 
   {
   startTimer(sensorTimer);
   for (unsigned ii=0; ii<sim.sensor_.size(); ++ii)
   {
      sim.sensor_[ii]->run(sim.time_, loop);
   }
   stopTimer(sensorTimer);
      
   startTimer(loopIOTimer);
   {
   if ( loop % sim.printRate_ == 0) 
   {
      int printIndex = sim.printIndex_; 
      if (printIndex >= 0)
      {
         printf("%8d %8.3f %12lld %21.15f %21.15f %21.15f\n",loop,sim.time_,sim.anatomy_.gid(printIndex),VmArray[printIndex],dVmR[printIndex],dVmD[printIndex]);  
         fprintf(sim.printFile_,"%8d %8.3f %12lld %21.15f %21.15f %21.15f\n",loop,sim.time_,sim.anatomy_.gid(printIndex),VmArray[printIndex],dVmR[printIndex],dVmD[printIndex]); 
      }
   }
   }

   if (sim.loop_ > 0 && sim.checkpointRate_ > 0 && sim.loop_ % sim.checkpointRate_ == 0)
      writeCheckpoint(sim, MPI_COMM_WORLD);

   if (!firstCall) 
   { 
      if (loop > 0 && loop % sim.snapshotRate_ == 0)
      {
         stringstream name;
         name << "snapshot."<<setfill('0')<<setw(12)<<loop;
         string fullname = name.str();
         if (myRank == 0) DirTestCreate(fullname.c_str());
         fullname += "/anatomy";
         writeCells(sim, fullname.c_str());
      }
   }
   
   }// critical section
   firstCall=0; 
   stopTimer(loopIOTimer);
}


void simulationLoop(Simulate& sim)
{
   vector<double> iStim(sim.anatomy_.nLocal(), 0.0);
   simulationProlog(sim);
   HaloExchange<double, AlignedAllocator<double> > voltageExchange(sim.sendMap_, (sim.commTable_));

   PotentialData& vdata=sim.vdata_;
   VectorDouble32& vmarray(vdata.VmArray_);
   VectorDouble32& dVmDiffusion(vdata.dVmDiffusion_);
   VectorDouble32& dVmReaction(vdata.dVmReaction_);

#if defined(SPI) && defined(TRACESPI)
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   cout << "Rank[" << myRank << "]: numOfNeighborToSend=" << sim.commTable_->_sendTask.size() << " numOfNeighborToRecv=" << sim.commTable_->_recvTask.size() << " numOfBytesToSend=" << sim.commTable_->_sendOffset[sim.commTable_->_sendTask.size()]*sizeof(double) << " numOfBytesToRecv=" << sim.commTable_->_recvOffset[sim.commTable_->_recvTask.size()]*sizeof(double) << endl;
#endif

   loopIO(sim,1);
   profileStart(simulationLoopTimer);
   while ( sim.loop_<sim.maxLoop_ )
   {
      int nLocal = sim.anatomy_.nLocal();
    
      startTimer(imbalanceTimer);
      voltageExchange.barrier();
      stopTimer(imbalanceTimer);
      
      startTimer(haloTimer);
      voltageExchange.fillSendBuffer(vmarray);
      voltageExchange.startComm();
      voltageExchange.wait();
      stopTimer(haloTimer);

      // DIFFUSION
      startTimer(diffusionCalcTimer);
      sim.diffusion_->updateLocalVoltage(&(vmarray[0]));
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
      sim.reaction_->calc(sim.dt_, vmarray, iStim, dVmReaction);
      stopTimer(reactionTimer);

      startTimer(integratorTimer);
      // no special BGQ integrator is this loop.  More bang for buck
      // from OMP threading.
      if (sim.checkRange_.on)
         sim.checkRanges(dVmReaction, dVmDiffusion);
      #pragma omp parallel for
      for (int ii=0; ii<nLocal; ++ii)
      {
         double dVm = dVmReaction[ii] + dVmDiffusion[ii];
         vmarray[ii] += sim.dt_*dVm;
      }

      sim.time_ += sim.dt_;
      ++sim.loop_;
      stopTimer(integratorTimer);
      
      if( sim.checkIO() )sim.bufferReactionData();
      
      loopIO(sim,0);
   }
   profileStop(simulationLoopTimer);
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
      timingBarrier = L2_BarrierWithSync_InitShared();

#ifdef PER_SQUAD_BARRIER
      {
	const int
	  nsq  = sim.reactionThreads_.nSquads(),
	  sqsz = sim.reactionThreads_.nThreads()/nsq /*sim.reactionThreads_.squadSize()*/;
	core_barrier = (L2_Barrier_t **) malloc(sizeof(L2_Barrier_t *) * nsq);
	for(int i = 0; i<nsq; i++)
	  core_barrier[i] = L2_BarrierWithSync_InitShared();
	integratorOffset_psb.resize(sim.reactionThreads_.nThreads() + 1);
	
	/*
	  This is an ugly import of the workBundle() function from
	   TT06Dev_Reaction.cc
	*/
	int workBundle(int index, int nItems, int nGroups , int mod, int& offset);
	integratorOffset_psb[0] = 0;
	for(int i = 0; i<nsq; i++) {
	  const int vlen = 4;
	  int ic,nc,nc_div,nc_mod;
	  nc = workBundle(i,sim.anatomy_.nLocal(),nsq,vlen,ic);
	  nc_div = nc / sqsz;
	  nc_mod = nc % sqsz;

	  if(integratorOffset_psb[i*sqsz] != ic) {
	    printf("%s:%d: Offset error, nsq=%d i=%d sqsz=%d nc=%d iOff[]=%d ic=%d\n",
		   __FILE__,__LINE__,nsq,i,sqsz,nc,integratorOffset_psb[i*sqsz],ic);
	    exit(1);
	  }

	  for(int j = 0; j<sqsz; j++)
	    integratorOffset_psb[i*sqsz + j + 1] = 
	      integratorOffset_psb[i*sqsz + j] + nc_div + (j < nc_mod);
	}
	  
      }
#endif

      int nLocal = sim.anatomy_.nLocal();
      const vector<int>& haloSendMap = voltageExchange.getSendMap();
      mkOffsets(integratorOffset, sim.anatomy_.nLocal(), sim.reactionThreads_);
      mkOffsets(fillSendBufOffset, haloSendMap.size(), sim.reactionThreads_);
   }

   ~SimLoopData()
   {
      free(reactionWaitOnNonGateBarrier);
      free(diffusionBarrier);
      free(reactionBarrier);
      free(timingBarrier);
   }
   
   
   L2_Barrier_t* haloBarrier;
   L2_Barrier_t* reactionBarrier;
   L2_Barrier_t* diffusionBarrier;
   L2_Barrier_t* reactionWaitOnNonGateBarrier;
   // The timing barrier syncs the top of the reaction loop with the
   // top of the diffusion loop.  This is done only to make it easier
   // to compare and understand timings.  This barrier can be removed
   // to slightly improve performance.
   L2_Barrier_t* timingBarrier;

#ifdef PER_SQUAD_BARRIER   
    L2_Barrier_t **core_barrier;
    vector<int> integratorOffset_psb;
#endif

   vector<int> integratorOffset;
   vector<pair<unsigned, unsigned> > sendBufIndirect;
   vector<int> fillSendBufOffset;
   HaloExchange<double, AlignedAllocator<double> > voltageExchange;
};
void diffusionLoopLag(Simulate& sim,
                   SimLoopData& loopData,
                   L2_BarrierHandle_t& reactionHandle,
                   L2_BarrierHandle_t& diffusionHandle)
{
   profileStart(diffusionLoopTimer);
   int tid = sim.diffusionThreads_.teamRank();
   L2_BarrierHandle_t haloBarrierHandle;
   L2_BarrierHandle_t timingHandle;
   L2_BarrierWithSync_InitInThread(loopData.haloBarrier, &haloBarrierHandle);
   L2_BarrierWithSync_InitInThread(loopData.timingBarrier, &timingHandle);
   int nTotalThreads = sim.reactionThreads_.nThreads() + sim.diffusionThreads_.nThreads();

   VectorDouble32& dVmDiffusion(sim.vdata_.dVmDiffusion_);

  // Need VmArray(t-h), Vm(t) and  remoteBuffer(t)    Note VmArray has local and remote and Vm just needs local. 

   while ( sim.loop_ < sim.maxLoop_ )
   {
      threadBarrier(timingBarrierTimer, loopData.timingBarrier, &timingHandle, nTotalThreads);
     
      // Halo Exchange
      if (tid == 0)
      {
         startTimer(haloLaunchTimer);
         loopData.voltageExchange.startComm();
         stopTimer(haloLaunchTimer);
         dVmDiffusion.assign(dVmDiffusion.size(), 0);
         for (unsigned ii=0; ii<sim.stimulus_.size(); ++ii)
            sim.stimulus_[ii]->stim(sim.time_, dVmDiffusion);
      }

      //stencil
      startTimer(diffusionCalcTimer);
      sim.diffusion_->calc(dVmDiffusion);
      stopTimer(diffusionCalcTimer);

      // announce that diffusion derivatives are ready.
      L2_BarrierWithSync_Arrive(loopData.diffusionBarrier, &diffusionHandle,
                                sim.diffusionThreads_.nThreads());
      L2_BarrierWithSync_Reset(loopData.diffusionBarrier,
                               &diffusionHandle,
                               sim.diffusionThreads_.nThreads());
      
      // Need a barrier for the completion of the halo exchange.
      if ( tid == 0) 
      {
         startTimer(haloWaitTimer);
         loopData.voltageExchange.wait();
         stopTimer(haloWaitTimer);
      }
      startTimer(diffusionL2BarrierHalo1Timer);
      L2_BarrierWithSync_Barrier(loopData.haloBarrier, &haloBarrierHandle, sim.diffusionThreads_.nThreads());
      stopTimer(diffusionL2BarrierHalo1Timer);
      
      // copy remote. Where is the timer?
      sim.diffusion_->updateRemoteVoltage(loopData.voltageExchange.getRecvBuf());

/********   With the timer barrier  this is not needed.  *******/

      // wait for reaction (integration) to finish
      startTimer(diffusionWaitTimer); 
      L2_BarrierWithSync_WaitAndReset(loopData.reactionBarrier,
                                      &reactionHandle,
                                      sim.reactionThreads_.nThreads());
      stopTimer(diffusionWaitTimer); 


   }
   profileStop(diffusionLoopTimer);
}

void diffusionLoop(Simulate& sim,
                   SimLoopData& loopData,
                   L2_BarrierHandle_t& reactionHandle,
                   L2_BarrierHandle_t& diffusionHandle)
{
   profileStart(diffusionLoopTimer);
   int tid = sim.diffusionThreads_.teamRank();
   L2_BarrierHandle_t haloBarrierHandle;
   L2_BarrierHandle_t timingHandle;
   L2_BarrierWithSync_InitInThread(loopData.haloBarrier, &haloBarrierHandle);
   L2_BarrierWithSync_InitInThread(loopData.timingBarrier, &timingHandle);
   int nTotalThreads = sim.reactionThreads_.nThreads() + sim.diffusionThreads_.nThreads();

   VectorDouble32& dVmDiffusion(sim.vdata_.dVmDiffusion_);


   uint64_t loopLocal = sim.loop_;
   while ( loopLocal < sim.maxLoop_ )
   {
     // Uncomment this for strong sync
      // if (tid == 0)
      // {
      //   startTimer(diffusionImbalanceTimer);
      //   loopData.voltageExchange.barrier();
      //   stopTimer(diffusionImbalanceTimer);
      // }
      threadBarrier(timingBarrierTimer, loopData.timingBarrier, &timingHandle, nTotalThreads);
     
      // Halo Exchange
      if (tid == 0)
      {
//         startTimer(diffusionImbalanceTimer);
//         loopData.voltageExchange.barrier();
//         stopTimer(diffusionImbalanceTimer);

         startTimer(haloTimer);

         startTimer(haloLaunchTimer);
         loopData.voltageExchange.startComm();
         stopTimer(haloLaunchTimer);

      // stimulus
         startTimer(stimulusTimer);
         dVmDiffusion.assign(dVmDiffusion.size(), 0);
         for (unsigned ii=0; ii<sim.stimulus_.size(); ++ii)
            sim.stimulus_[ii]->stim(sim.time_, dVmDiffusion);
         stopTimer(stimulusTimer);

         startTimer(haloWaitTimer);
         loopData.voltageExchange.wait();
         stopTimer(haloWaitTimer);

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
      sim.diffusion_->calc(dVmDiffusion);
      stopTimer(diffusionCalcTimer);
      
      //barrier
      startTimer(diffusionL2BarrierHalo2Timer);
      L2_BarrierWithSync_Barrier(loopData.haloBarrier, &haloBarrierHandle,
                                 sim.diffusionThreads_.nThreads());
      stopTimer(diffusionL2BarrierHalo2Timer);


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
      ++loopLocal;
//      startTimer(dummyTimer);
//      stopTimer(dummyTimer);

   }
   profileStop(diffusionLoopTimer);
}
void reactionLoop(Simulate& sim, SimLoopData& loopData, L2_BarrierHandle_t& reactionHandle, L2_BarrierHandle_t& diffusionHandle, 
                  void (*integrate)(const int, const int, const double, double*, double*, unsigned*, double*, double*k, double*, double))
{
   profileStart(reactionLoopTimer);
   int tid = sim.reactionThreads_.teamRank();
   L2_BarrierHandle_t reactionWaitOnNonGateHandle;
   L2_BarrierHandle_t timingHandle;
   L2_BarrierWithSync_InitInThread(loopData.reactionWaitOnNonGateBarrier, &reactionWaitOnNonGateHandle);
   L2_BarrierWithSync_InitInThread(loopData.timingBarrier, &timingHandle);
   VectorDouble32& VmArray(sim.vdata_.VmArray_);
   VectorDouble32& dVmReaction(sim.vdata_.dVmReaction_);
   VectorDouble32& dVmDiffusion(sim.vdata_.dVmDiffusion_);
   int nTotalThreads = sim.reactionThreads_.nThreads() + sim.diffusionThreads_.nThreads();

#ifdef PER_SQUAD_BARRIER
   const int sqsz = sim.reactionThreads_.squadSize();
   const int b_id = sim.reactionThreads_.rankInfo().coreRank_;
   L2_BarrierHandle_t core_barrier_h;
   L2_Barrier_t *cb_ptr = loopData.core_barrier[b_id];
   L2_BarrierWithSync_InitInThread(cb_ptr, &core_barrier_h);
#endif


   if (tid == 0) 
   {
     loopIO(sim,1);
   }
   uint64_t loopLocal = sim.loop_;
   while ( loopLocal < sim.maxLoop_ )
   {
      threadBarrier(timingBarrierTimer, loopData.timingBarrier, &timingHandle, nTotalThreads);
           
      startTimer(reactionTimer);
      sim.reaction_->updateNonGate(sim.dt_, VmArray, dVmReaction);

      startTimer(GateNonGateTimer); 
#ifdef PER_SQUAD_BARRIER
      {
	L2_Barrier_nosync_Arrive(      cb_ptr,
				       &core_barrier_h,
				       sqsz);
	L2_Barrier_nosync_WaitAndReset(cb_ptr,
				       &core_barrier_h,
				       sqsz);
      }
#else
      L2_BarrierWithSync_Barrier(loopData.reactionWaitOnNonGateBarrier,
				 &reactionWaitOnNonGateHandle,
				 sim.reactionThreads_.nThreads());
#endif
      stopTimer(GateNonGateTimer); 
      sim.reaction_->updateGate(sim.dt_, VmArray);
      stopTimer(reactionTimer);

      startTimer(reactionWaitTimer);
      /*
	    The following barrier makes sure the integrator does not write to data still being read by some
	    threads still in updateGate(). With the modified loop order in integratorOffset_psb[], only per
	    squad (per core) barriers are needed. With the original loop order in integratorOffset[], a full
	    barrier over the reaction cores is necessary.  
      */
#ifdef PER_SQUAD_BARRIER
      {
	      L2_Barrier_nosync_Arrive(      cb_ptr, &core_barrier_h, sqsz);
	      L2_Barrier_nosync_WaitAndReset(cb_ptr, &core_barrier_h, sqsz);
      }
#else
      L2_BarrierWithSync_Barrier(loopData.reactionWaitOnNonGateBarrier,
				 &reactionWaitOnNonGateHandle,
				 sim.reactionThreads_.nThreads());
#endif

      L2_BarrierWithSync_WaitAndReset(loopData.diffusionBarrier, &diffusionHandle, sim.diffusionThreads_.nThreads());
      stopTimer(reactionWaitTimer);

      startTimer(dummyTimer);
      stopTimer(dummyTimer);

      startTimer(integratorTimer);
#ifdef PER_SQUAD_BARRIER
      const int begin=loopData.integratorOffset_psb[tid];
      const int end  =loopData.integratorOffset_psb[tid+1]
#else
      const int begin=loopData.integratorOffset[tid];
      const int end  =loopData.integratorOffset[tid+1];
#endif

      integrate(begin, end, sim.dt_, &dVmReaction[0], &dVmDiffusion[0],
		    sim.diffusion_->blockIndex(),
		    sim.diffusion_->dVmBlock(),
		    sim.diffusion_->VmBlock(),
		    &VmArray[0],
		    sim.diffusion_->diffusionScale());
      stopTimer(integratorTimer);

      startTimer(reactionMiscTimer); 
      if (sim.checkRange_.on) sim.checkRanges(begin, end, dVmReaction, dVmDiffusion);
      ++loopLocal;
      if (tid == 0)
      {
         sim.time_ += sim.dt_;
         ++sim.loop_;
      }
      
      if ( sim.checkIO(loopLocal) ) // every thread should return the same result
      {
         // jlf: needed to have correct dVm in sensors.
         // Already calculated in integrateLoop, so maybe we 
         // should save it there and reuse it here...

         const unsigned* const blockIndex = sim.diffusion_->blockIndex();
         const double* const dVdMatrix = sim.diffusion_->dVmBlock();
         const double diffusionScale = sim.diffusion_->diffusionScale();
         for (unsigned ii=begin; ii<end; ++ii)
         {
            const int index = blockIndex[ii];
            // add diffusion to stimulus (already stored in dVmDiffusion)
            dVmDiffusion[ii] += diffusionScale*dVdMatrix[index];
         }
            
         sim.bufferReactionData(begin, end);
      }
         
      L2_BarrierWithSync_Barrier(loopData.reactionWaitOnNonGateBarrier,
                                 &reactionWaitOnNonGateHandle,
                                 sim.reactionThreads_.nThreads());
      if (tid == 0) loopIO(sim,0);
      stopTimer(reactionMiscTimer); 
         
      startTimer(haloMove2BufTimer);
      {
	const vector<int>& haloSendMap = loopData.voltageExchange.getSendMap();
         unsigned begin = loopData.fillSendBufOffset[tid];
         unsigned end   = loopData.fillSendBufOffset[tid+1];
         double* sendBuf = loopData.voltageExchange.getSendBuf();
         for (unsigned ii=begin; ii<end; ++ii)
            sendBuf[ii] = VmArray[haloSendMap[ii]];
      }
      stopTimer(haloMove2BufTimer);
      
      startTimer(reactionL2ArriveTimer);
      L2_BarrierWithSync_Arrive(loopData.reactionBarrier, &reactionHandle, sim.reactionThreads_.nThreads());
      stopTimer(reactionL2ArriveTimer);
      
      startTimer(reactionL2ResetTimer);
      L2_BarrierWithSync_Reset(loopData.reactionBarrier, &reactionHandle, sim.reactionThreads_.nThreads());
      stopTimer(reactionL2ResetTimer);
      //#pragma omp barrier 
   }
   profileStop(reactionLoopTimer);
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
   VectorDouble32& VmArray(sim.vdata_.VmArray_);

   timestampBarrier("Entering threaded region", MPI_COMM_WORLD);

   #pragma omp parallel
   {
      int ompTid = omp_get_thread_num();
      
      L2_BarrierHandle_t reactionHandle;
      L2_BarrierHandle_t diffusionHandle;
      L2_BarrierWithSync_InitInThread(loopData.reactionBarrier, &reactionHandle);
      L2_BarrierWithSync_InitInThread(loopData.diffusionBarrier, &diffusionHandle);
      
      // setup matrix voltages for first timestep.
      if (sim.reactionThreads_.teamRank() >= 0) 
         sim.diffusion_->updateLocalVoltage(&VmArray[0]);
      if (sim.reactionThreads_.teamRank() == 0) 
         loopData.voltageExchange.fillSendBuffer(VmArray);
      #pragma omp barrier
      profileStart(simulationLoopTimer);
      if ( sim.diffusionThreads_.teamRank() >= 0) 
      {
         diffusionLoop(sim, loopData, reactionHandle, diffusionHandle);
      }
      if ( sim.reactionThreads_.teamRank() >= 0) 
      {
          reactionLoop(sim, loopData, reactionHandle, diffusionHandle,integrateLoop);
      } 
      profileStop(simulationLoopTimer);
   }
}
void simulationLoopLag(Simulate& sim)
{
   SimLoopData loopData(sim);

   simulationProlog(sim);

#if defined(SPI) && defined(TRACESPI)
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   cout << "Rank[" << myRank << "]: numOfNeighborToSend=" << sim.commTable_->_sendTask.size() << " numOfNeighborToRecv=" << sim.commTable_->_recvTask.size() << " numOfBytesToSend=" << sim.commTable_->_sendOffset[sim.commTable_->_sendTask.size()]*sizeof(double) << " numOfBytesToRecv=" << sim.commTable_->_recvOffset[sim.commTable_->_recvTask.size()]*sizeof(double) << endl;
#endif
   VectorDouble32& VmArray(sim.vdata_.VmArray_);

   timestampBarrier("Entering threaded region", MPI_COMM_WORLD);

   #pragma omp parallel
   {
      int ompTid = omp_get_thread_num();
      
      L2_BarrierHandle_t reactionHandle;
      L2_BarrierHandle_t diffusionHandle;
      L2_BarrierWithSync_InitInThread(loopData.reactionBarrier, &reactionHandle);
      L2_BarrierWithSync_InitInThread(loopData.diffusionBarrier, &diffusionHandle);
      
      // setup matrix voltages for first timestep.
      if (sim.reactionThreads_.teamRank() >= 0) 
      {
         sim.diffusion_->updateLocalVoltage(&VmArray[0]);      
         loopData.voltageExchange.fillSendBuffer(VmArray);     // Check which tasks need to call fillSendBuffer JNG
      }
      #pragma omp barrier
      if ( sim.diffusionThreads_.teamRank() == 0) // Boot strap step.  Don't know Vm(t-h) use Vm(t) and 
      {
           loopData.voltageExchange.startComm();
           loopData.voltageExchange.wait(); 
           sim.diffusion_->updateRemoteVoltage(loopData.voltageExchange.getRecvBuf());
      }
      #pragma omp barrier
      profileStart(simulationLoopTimer);
      if ( sim.diffusionThreads_.teamRank() >= 0) 
      {
         diffusionLoopLag(sim, loopData, reactionHandle, diffusionHandle);
      }
      if ( sim.reactionThreads_.teamRank() >= 0) 
      {
          reactionLoop(sim, loopData, reactionHandle, diffusionHandle,integrateLoopLag);   //Need to change Integration Routine to update VmArray with Vm(t) and not Vm(t+h); 
      } 
      profileStop(simulationLoopTimer);
   }
}

void simulationLoopAllSkate(Simulate& sim)
{
   L2_Barrier_t* barrierPtr = L2_BarrierWithSync_InitShared();
   simulationProlog(sim);
   HaloExchange<double, AlignedAllocator<double> >
      voltageExchange(sim.sendMap_, (sim.commTable_));
   loopIO(sim, 1);

   // define useful aliases
   VectorDouble32& VmArray(sim.vdata_.VmArray_);
   VectorDouble32& dVmReaction(sim.vdata_.dVmReaction_);
   VectorDouble32& dVmDiffusion(sim.vdata_.dVmDiffusion_);
   const ThreadTeam& threadInfo = sim.reactionThreads_;
   const unsigned nLocal = sim.anatomy_.nLocal();
   const vector<int>& haloSendMap = voltageExchange.getSendMap();
   double* haloSendBuf = voltageExchange.getSendBuf();
   int nThreads = threadInfo.nThreads();
   
   
   vector<int> fillSendBufOffset(nThreads+1);
   vector<int> zeroDiffusionOffset(nThreads+1);
   vector<int> integratorOffset(nThreads+1);
   mkOffsets(fillSendBufOffset, haloSendMap.size(), threadInfo);
   mkOffsets(zeroDiffusionOffset, nLocal, threadInfo);
   mkOffsets(integratorOffset, nLocal, threadInfo);

   timestampBarrier("Entering threaded region", MPI_COMM_WORLD);
   #pragma omp parallel
   {
      int tid = threadInfo.teamRank(); 
      int fsboBegin = fillSendBufOffset[tid];
      int fsboEnd =   fillSendBufOffset[tid+1];
      int zdoBegin = zeroDiffusionOffset[tid];
      int zdoEnd =   zeroDiffusionOffset[tid+1];
      int ioBegin = integratorOffset[tid];
      int ioEnd =   integratorOffset[tid+1];
      
      L2_BarrierHandle_t barrierHandle;
      L2_BarrierWithSync_InitInThread(barrierPtr, &barrierHandle);

      startTimer(haloMove2BufTimer);
      for (int ii = fsboBegin; ii<fsboEnd; ++ii)
         haloSendBuf[ii] = VmArray[haloSendMap[ii]];
      stopTimer(haloMove2BufTimer);
      sim.diffusion_->updateLocalVoltage(&VmArray[0]);

     #pragma omp barrier
      
      profileStart(simulationLoopTimer);
      while (sim.loop_ < sim.maxLoop_)
      {
         // This is a big barrier to catch all imbalance if we want.
//          threadBarrier(threadImbalanceTimer, barrierPtr, &barrierHandle, nThreads)
//          startTimer(nodeImbalanceTimer);
//          voltageExchange.barrier();
//          stopTimer(NodeImbalanceTimer);

         if (tid == 0)
         {
            startTimer(haloLaunchTimer);
            voltageExchange.startComm();
            stopTimer(haloLaunchTimer);
         }
         
         sim.reaction_->updateNonGate(sim.dt_, VmArray, dVmReaction);
         
         threadBarrier(barrier1Timer, barrierPtr, &barrierHandle, nThreads);

         sim.reaction_->updateGate(sim.dt_, VmArray);

         if (tid == 0)
         {
            startTimer(haloWaitTimer);
            voltageExchange.wait();
            stopTimer(haloWaitTimer);
         }
         
         threadBarrier(barrier1Timer, barrierPtr, &barrierHandle, nThreads);

         sim.diffusion_->updateRemoteVoltage(voltageExchange.getRecvBuf());

         startTimer(initializeDVmDTimer);
         for (int ii=zdoBegin; ii<zdoEnd; ++ii)
            dVmDiffusion[ii] = 0;
         stopTimer(initializeDVmDTimer);
         
         threadBarrier(barrier2Timer, barrierPtr, &barrierHandle, nThreads);

         if (tid == 0)
         {
            startTimer(stimulusTimer);
            for (unsigned ii=0; ii<sim.stimulus_.size(); ++ii)
               sim.stimulus_[ii]->stim(sim.time_, dVmDiffusion);
            stopTimer(stimulusTimer);
         }
         
         threadBarrier(barrier3Timer, barrierPtr, &barrierHandle, nThreads);
         
         startTimer(diffusionCalcTimer);
         sim.diffusion_->calc(dVmDiffusion);
         stopTimer(diffusionCalcTimer);

         threadBarrier(barrier4Timer, barrierPtr, &barrierHandle, nThreads);

         startTimer(integratorTimer);
         integrateLoop(ioBegin, ioEnd, sim.dt_, &dVmReaction[0], &dVmDiffusion[0],
                       sim.diffusion_->blockIndex(),
                       sim.diffusion_->dVmBlock(),
                       sim.diffusion_->VmBlock(),
                       &VmArray[0],
                       sim.diffusion_->diffusionScale());
         stopTimer(integratorTimer);

         if (sim.checkRange_.on)
         {
            startTimer(rangeCheckTimer); 
            sim.checkRanges(ioBegin, ioEnd, dVmReaction, dVmDiffusion);
            stopTimer(rangeCheckTimer); 
         }

         if (tid == 0)
         {
            sim.time_ += sim.dt_;
            ++sim.loop_;
         }
         
         if ( tid == 0 && sim.loop_ % sim.printRate_ == 0) 
         {
            startTimer(printDataTimer);
            int pi = sim.printIndex_; 
            if (pi >= 0)
            {
               const unsigned* const blockIndex = sim.diffusion_->blockIndex();
               const double* const dVdMatrix = sim.diffusion_->dVmBlock();
               const double scale = sim.diffusion_->diffusionScale();
               const int index = blockIndex[pi];
               double dVd = dVmDiffusion[pi] + scale*dVdMatrix[index];
               printf("%8d %8.3f %12lld %21.15f %21.15f %21.15f\n",
                      sim.loop_, sim.time_, sim.anatomy_.gid(pi),
                      VmArray[pi], dVmReaction[pi], dVd);  
               fprintf(sim.printFile_,
                       "%8d %8.3f %12lld %21.15f %21.15f %21.15f\n",
                       sim.loop_, sim.time_, sim.anatomy_.gid(pi),
                       VmArray[pi], dVmReaction[pi], dVd);  
            }
            stopTimer(printDataTimer);
         }
         
//         sensors;
//         checkpoint;

         startTimer(haloMove2BufTimer);
         for (unsigned ii = fsboBegin; ii<fsboEnd; ++ii)
            haloSendBuf[ii] = VmArray[haloSendMap[ii]];
         stopTimer(haloMove2BufTimer);

         threadBarrier(barrier5Timer, barrierPtr, &barrierHandle, nThreads);
      }
      profileStop(simulationLoopTimer);
   }
   //   free(barrierPtr);
}
