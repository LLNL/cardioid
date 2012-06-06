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
   void mkSendBufInfo(vector<unsigned>& sendBufOffset,
                      vector<pair<unsigned, unsigned> >& sendBufIndirect,
                      const HaloExchange<double>& halo,
                      const vector<int>& integratorOffset,
                      const ThreadTeam& threadInfo)
   {
      int nThreads = threadInfo.nThreads();
      sendBufOffset.resize(nThreads+1, 0);
      const vector<int>& sendMap = halo.getSendMap();
      sendBufIndirect.reserve(sendMap.size());
      
      for (unsigned ii=0; ii<sendMap.size(); ++ii)
         sendBufIndirect.push_back(make_pair(sendMap[ii], ii));
      sort(sendBufIndirect.begin(), sendBufIndirect.end());
      
      for (unsigned ii=0; ii<nThreads; ++ii)
      {
         unsigned target = integratorOffset[ii+1];
         for (unsigned jj=sendBufOffset[ii]; jj<sendBufIndirect.size(); ++jj)
         {
            if (sendBufIndirect[jj].first >= target)
            {
               sendBufOffset[ii+1] = jj;
               break;
            }
         }
         if (sendBufOffset[ii+1] == 0) // didn't find an end point
            sendBufOffset[ii+1] = sendBufIndirect.size();
      }
   }
}


void simulationProlog(Simulate& sim)
{
   // initialize membrane voltage with default value from the reaction model. 
   sim.vdata_.setup(sim.anatomy_);
   //sim.VmArray_.resize(sim.anatomy_.size());
   
   std::vector<double>& VmArray(*sim.vdata_.VmArray_);
   sim.reaction_->initializeMembraneVoltage(VmArray);
  
   for (unsigned ii=sim.anatomy_.nLocal(); ii<sim.anatomy_.size(); ++ii)
      VmArray[ii] = 0;

   // Load state file, assign corresponding values to membrane voltage and cell model
   for (unsigned ii=0; ii<sim.stateFilename_.size(); ++ii)
      readCheckpoint(sim.stateFilename_[ii], sim, MPI_COMM_WORLD);

}


void loopIO(const Simulate& sim)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
      
   static int firstCall=1; 
   int loop = sim.loop_; 

   std::vector<double>& VmArray(*sim.vdata_.VmArray_);
   std::vector<double>& dVmR(*sim.vdata_.dVmReaction_);
   std::vector<double>& dVmD(*sim.vdata_.dVmDiffusion_);

   // SENSORS
   startTimer(sensorTimer);
   for (unsigned ii=0; ii<sim.sensor_.size(); ++ii)
   {
       sim.sensor_[ii]->run(sim.time_, loop);
   }
   stopTimer(sensorTimer);
      
   startTimer(loopIOTimer);
   if ( (loop % sim.printRate_ == 0) && myRank == 0)
   {
      static FILE *file; 
      if (firstCall) 
      {
      file=fopen("data","a"); 
       printf("#   Loop     Time         Vm(t)        dVm_r(t-h)      dVm_d(t-h)\n");
       fprintf(file,"#   Loop     Time         Vm(t)        dVm_r(t-h)      dVm_d(t-h)\n");
      }
      printf("%8d %8.3f %21.15f %21.15f %21.15f\n",loop,sim.time_,VmArray[0],dVmR[0],dVmD[0]);  fflush(stdout); 
      fprintf(file,"%8d %8.3f %21.15f %21.15f %21.15f\n",loop,sim.time_,VmArray[0],dVmR[0],dVmD[0]); fflush(file); 
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
   firstCall=0; 
   stopTimer(loopIOTimer);
}


void simulationLoop(Simulate& sim)
{
   vector<double> iStim(sim.anatomy_.nLocal(), 0.0);
   simulationProlog(sim);
   HaloExchange<double> voltageExchange(sim.sendMap_, (sim.commTable_));

   PotentialData& vdata=sim.vdata_;
   std::vector<double>& vmarray(*vdata.VmArray_);
   std::vector<double>& dVmDiffusion(*vdata.dVmDiffusion_);
   std::vector<double>& dVmReaction(*vdata.dVmReaction_);

#if defined(SPI) && defined(TRACESPI)
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   cout << "Rank[" << myRank << "]: numOfNeighborToSend=" << sim.commTable_->_sendTask.size() << " numOfNeighborToRecv=" << sim.commTable_->_recvTask.size() << " numOfBytesToSend=" << sim.commTable_->_sendOffset[sim.commTable_->_sendTask.size()]*sizeof(double) << " numOfBytesToRecv=" << sim.commTable_->_recvOffset[sim.commTable_->_recvTask.size()]*sizeof(double) << endl;
#endif

   loopIO(sim);
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
      if (sim.checkRanges_)
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
      loopIO(sim);
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
      dVmReaction=new vector<double>(nLocal, 0.0);
      mkOffsets(integratorOffset, sim.anatomy_.nLocal(), sim.reactionThreads_);
      mkSendBufInfo(sendBufOffset, sendBufIndirect, voltageExchange,
                    integratorOffset, sim.reactionThreads_);
   }

   ~SimLoopData()
   {
      free(reactionWaitOnNonGateBarrier);
      free(diffusionBarrier);
      free(reactionBarrier);
      delete dVmReaction;
   }
   
   
   L2_Barrier_t* haloBarrier;
   L2_Barrier_t* reactionBarrier;
   L2_Barrier_t* diffusionBarrier;
   L2_Barrier_t* reactionWaitOnNonGateBarrier;

#ifdef PER_SQUAD_BARRIER   
    L2_Barrier_t **core_barrier;
    vector<int> integratorOffset_psb;
#endif

   vector<double>* dVmReaction; 
   vector<int> integratorOffset;
   vector<pair<unsigned, unsigned> > sendBufIndirect;
   vector<unsigned> sendBufOffset;
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

   std::vector<double>& dVmDiffusion(*sim.vdata_.dVmDiffusion_);
   std::vector<double>& dVmReaction(*sim.vdata_.dVmReaction_);

   if (tid == 0)
      loopIO(sim);
   while ( sim.loop_ < sim.maxLoop_ )
   {
      int nLocal = sim.anatomy_.nLocal();
      
      // Halo Exchange
      if (tid == 0)
      {
         startTimer(diffusionImbalanceTimer);
//         loopData.voltageExchange.barrier();
         stopTimer(diffusionImbalanceTimer);

         startTimer(haloTimer);

         startTimer(haloTimerExecute);
         loopData.voltageExchange.startComm();
         stopTimer(haloTimerExecute);

      // stimulus
         startTimer(stimulusTimer);
         dVmDiffusion.assign(dVmDiffusion.size(), 0);
         for (unsigned ii=0; ii<sim.stimulus_.size(); ++ii)
            sim.stimulus_[ii]->stim(sim.time_, dVmDiffusion);
         stopTimer(stimulusTimer);

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

      startTimer(dummyTimer);
      stopTimer(dummyTimer);

      
      //#pragma omp barrier
      if (tid == 0 && sim.loop_ % sim.printRate_ == 0)
      {
         startTimer(diffusiondVmRCopyTimer);
#if 0
         //jlf: code to swap pointers and avoid copy. 
         //Not working... may need a barrier to sync with reaction loop?
         loopData.dVmReaction = sim.vdata_.swapdVmReaction(loopData.dVmReaction);
         dVmReaction=(*sim.vdata_.dVmReaction_); // update reference
#else
         dVmReaction = *loopData.dVmReaction; 
#endif
         stopTimer(diffusiondVmRCopyTimer);

         loopIO(sim);
      }
   }
   profileStop(diffusionLoopTimer);
}

void reactionLoop(Simulate& sim, SimLoopData& loopData, L2_BarrierHandle_t& reactionHandle, L2_BarrierHandle_t& diffusionHandle)
{
   profileStart(reactionLoopTimer);
   int tid = sim.reactionThreads_.teamRank();
   L2_BarrierHandle_t reactionWaitOnNonGateHandle;
   L2_BarrierWithSync_InitInThread(loopData.reactionWaitOnNonGateBarrier, &reactionWaitOnNonGateHandle);
   std::vector<double>& VmArray(*sim.vdata_.VmArray_);
   std::vector<double>& dVmDiffusion(*sim.vdata_.dVmDiffusion_);

#ifdef PER_SQUAD_BARRIER
   const int sqsz = sim.reactionThreads_.squadSize();
   const int b_id = sim.reactionThreads_.rankInfo().coreRank_;
   L2_BarrierHandle_t core_barrier_h;
   L2_Barrier_t *cb_ptr = loopData.core_barrier[b_id];
   L2_BarrierWithSync_InitInThread(cb_ptr, &core_barrier_h);
#endif

   while ( sim.loop_ < sim.maxLoop_ )
   {
      int nLocal = sim.anatomy_.nLocal();
      
      startTimer(reactionTimer);
      sim.reaction_->updateNonGate(sim.dt_, VmArray, *loopData.dVmReaction);
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


      /*
	The following barrier makes sure the integrator
	does not write to data still beaing read by some
	threads still in updateGate(). With the modified
	loop order in integratorOffset_psb[], only per
	squad (per core) barriers are needed. With the
	original loop order in integratorOffset[], a full
	barrier over the reaction cores is necessary.
      */
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


      startTimer(dummyTimer);
      stopTimer(dummyTimer);

      startTimer(reactionWaitTimer);
      L2_BarrierWithSync_WaitAndReset(loopData.diffusionBarrier, &diffusionHandle, sim.diffusionThreads_.nThreads());
      stopTimer(reactionWaitTimer);


//       copy the diffusion derivative from matrix to array form.
//       This is a dinosaur track.  The code that actually does this
//       work is now inside the integrator.      
//       startTimer(FGR_Matrix2ArrayTimer);
//       unsigned begin = loopData.integratorOffset[tid];
//       unsigned end   = loopData.integratorOffset[tid+1];
//       unsigned* blockIndex = sim.diffusion_->blockIndex();
//       double* dVdMatrix = sim.diffusion_->dVmBlock();
//       double diffusionScale = sim.diffusion_->diffusionScale();
//       for (unsigned ii=begin; ii<end; ++ii)
//       {
//          int index = blockIndex[ii];
//          loopData.dVmDiffusion[ii] += diffusionScale*dVdMatrix[index];
//       }
//       stopTimer(FGR_Matrix2ArrayTimer);

      
      startTimer(integratorTimer);
#ifdef PER_SQUAD_BARRIER
      integrateLoop(loopData.integratorOffset_psb[tid],
		    loopData.integratorOffset_psb[tid+1],
		    sim.dt_,
		    &(*loopData.dVmReaction)[0],
		    &dVmDiffusion[0],
		    sim.diffusion_->blockIndex(),
		    sim.diffusion_->blockIndexB(),
		    sim.diffusion_->dVmBlock(),
		    sim.diffusion_->VmBlock(),
		    &VmArray[0],
		    sim.diffusion_->diffusionScale());
#else
      integrateLoop(loopData.integratorOffset[tid],
		    loopData.integratorOffset[tid+1],
		    sim.dt_,
		    &(*loopData.dVmReaction)[0],
		    &dVmDiffusion[0],
		    sim.diffusion_->blockIndex(),
		    sim.diffusion_->blockIndexB(),
		    sim.diffusion_->dVmBlock(),
		    sim.diffusion_->VmBlock(),
		    &VmArray[0],
		    sim.diffusion_->diffusionScale());
#endif
      
      if (tid == 0)
      {
         if (sim.checkRanges_)
            sim.checkRanges(*loopData.dVmReaction, dVmDiffusion);
         
         sim.time_ += sim.dt_;
         ++sim.loop_;
      }

      stopTimer(integratorTimer);
//       sim.diffusion_->updateLocalVoltage(&sim.VmArray_[0]);

      startTimer(haloMove2BufTimer);
      {
         unsigned begin = loopData.sendBufOffset[tid];
         unsigned end   = loopData.sendBufOffset[tid+1];
         double* sendBuf = loopData.voltageExchange.getSendBuf();
         for (unsigned ii=begin; ii<end; ++ii)
         {
            unsigned vIndex = loopData.sendBufIndirect[ii].first;
            unsigned sIndex = loopData.sendBufIndirect[ii].second;
            sendBuf[sIndex] = VmArray[vIndex];
         }
      }
      stopTimer(haloMove2BufTimer);
      
      startTimer(reactionL2ArriveTimer);
      L2_BarrierWithSync_Arrive(loopData.reactionBarrier, &reactionHandle, sim.reactionThreads_.nThreads());
      stopTimer(reactionL2ArriveTimer);

      startTimer(reactionL2ResetTimer);
      L2_BarrierWithSync_WaitAndReset(loopData.reactionBarrier, &reactionHandle, sim.reactionThreads_.nThreads());
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
   std::vector<double>& VmArray(*sim.vdata_.VmArray_);

   #pragma omp parallel
   {
      int ompTid = omp_get_thread_num();
      
      L2_BarrierHandle_t reactionHandle;
      L2_BarrierHandle_t diffusionHandle;
      L2_BarrierWithSync_InitInThread(loopData.reactionBarrier, &reactionHandle);
      L2_BarrierWithSync_InitInThread(loopData.diffusionBarrier, &diffusionHandle);
      
      // setup matrix voltages for first timestep.
      sim.diffusion_->updateLocalVoltage(&VmArray[0]);
      loopData.voltageExchange.fillSendBuffer(VmArray);
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
