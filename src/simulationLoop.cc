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
#include "PerformanceTimers.hh"
#include "BucketOfBits.hh"
#include "stateLoader.hh"
#include "fastBarrier.hh"

using namespace std;


static const int nDiffuse = 1;
static int nReact =0; 
static L2_Barrier_t reactionBarrier;
static L2_Barrier_t diffusionBarrier;
   
void simulationProlog(Simulate& sim, vector<double> dVmReaction, vector < double> dVmDiffusion, vector < double> dVmExternal,int myRank)
{
  // initialize membrane voltage with a default value from the reaction
  // model. 
  sim.VmArray_.resize(sim.anatomy_.size());
  sim.reaction_->initializeMembraneVoltage(sim.VmArray_);
  
  for (unsigned ii=sim.anatomy_.nLocal(); ii<sim.anatomy_.size(); ++ii)
     sim.VmArray_[ii] = 0;

//ddt  for (unsigned ii=0; ii<sim.anatomy_.nLocal(); ++ii)
//ddt  cout << sim.anatomy_.gid(ii) <<endl;
  
  // Load state file, assign corresponding values to membrane voltage
  // and cell model
  if (!sim.stateFilename_.empty())
  {
     BucketOfBits stateData = loadAndDistributeState(sim.stateFilename_, sim.anatomy_);
     assert(stateData.nRecords() == sim.anatomy_.nLocal());
     sim.reaction_->loadState(stateData);
     unsigned vmIndex = stateData.getIndex("Vm");
     if (vmIndex != stateData.nFields())
        for (unsigned ii=0; ii<stateData.nRecords(); ++ii)
           stateData.getRecord(ii).getValue(vmIndex, sim.VmArray_[ii]);
  }
  if ( myRank == 0)
  {
    cout << "    Loop     Time           Vm        dVm_r        dVm_d        dVm_e" <<endl;
    cout << setw(8) << sim.loop_ <<" "
         << setw(8) << sim.time_ <<" "
         << setw(12) << sim.VmArray_[0] << " "
         << setw(12) << dVmReaction[0] << " "
         << setw(12) << dVmDiffusion[0] << " "
         << setw(12) << dVmExternal[0]  << endl;
  }
  // sensors:  print out initial values
  for (unsigned ii=0; ii<sim.sensor_.size(); ++ii)
  {
     sim.sensor_[ii]->print(sim.time_, sim.loop_, sim.VmArray_, dVmReaction, dVmDiffusion, dVmExternal);
  }
}
void loopIO(Simulate& sim, vector<double> dVmReaction, vector < double> dVmDiffusion, vector < double> dVmExternal,int myRank)
{
    
    static TimerHandle loopIOHandle = profileGetHandle("LoopIO");
    profileStart(loopIOHandle);

    for (unsigned ii=0; ii<sim.sensor_.size(); ++ii)
    {
      if (sim.loop_ % sim.sensor_[ii]->printRate() == 0)
         sim.sensor_[ii]->print(sim.time_, sim.loop_, sim.VmArray_, dVmReaction, dVmDiffusion, dVmExternal);
    }
    
    if ( (sim.loop_ % sim.printRate_ == 0) && myRank == 0)
    {
      cout << setw(8) << sim.loop_ <<" "
           << setw(8) << sim.time_ <<" "
           << setw(12) << sim.VmArray_[0] << " "
           << setw(12) << dVmReaction[0] << " "
           << setw(12) << dVmDiffusion[0] << " "
           << setw(12) << dVmExternal[0]  << " "
           << setw(12) << sim.anatomy_.gid(0)  << endl;
    }
    
    if (sim.loop_ % sim.snapshotRate_ == 0)
    {
      stringstream name;
      name << "snapshot."<<setfill('0')<<setw(8)<<sim.loop_;
      string fullname = name.str();
      if (myRank == 0)
         DirTestCreate(fullname.c_str());
      fullname += "/anatomy";
      writeCells(sim, fullname.c_str());
    }
    profileStop(loopIOHandle);
}

void simulationLoop(Simulate& sim)
{
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  vector<double> dVmDiffusion(sim.anatomy_.nLocal(), 0.0);
  vector<double> dVmReaction(sim.anatomy_.nLocal(), 0.0);
  vector<double> dVmExternal(sim.anatomy_.nLocal(), 0.0);
  vector<double> iStim(sim.anatomy_.nLocal(), 0.0);
  simulationProlog(sim, dVmReaction,dVmDiffusion,dVmExternal,myRank);
#ifdef SPI   
  spi_HaloExchange<double> voltageExchange(sim.sendMap_, (sim.commTable_));
#else
  mpi_HaloExchange<double> voltageExchange(sim.sendMap_, (sim.commTable_));
#endif
  while ( sim.loop_<=sim.maxLoop_ )
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

    // SENSORS
    static TimerHandle sensorHandle = profileGetHandle("Sensors");
    profileStart(sensorHandle);
    for (unsigned ii=0; ii<sim.sensor_.size(); ++ii)
    {
      if (sim.loop_ % sim.sensor_[ii]->evalRate() == 0)
         sim.sensor_[ii]->eval(sim.time_, sim.loop_, sim.VmArray_, dVmReaction, dVmDiffusion, dVmExternal);
    }
    profileStop(sensorHandle);
    loopIO(sim, dVmReaction, dVmDiffusion, dVmExternal, myRank);
  }
}
void diffusionLoop(Simulate& sim, vector<double>& dVmReaction, int myRank, L2_BarrierHandle_t& reactionHandle, L2_BarrierHandle_t& diffusionHandle)
{
  int ompTid = omp_get_thread_num();
  printf("diffusionLoop: ompTid=%d\n",ompTid); 
  vector<double> dVm(sim.anatomy_.nLocal(), 0.0);
  vector<double> iStim(sim.anatomy_.nLocal(), 0.0);
#ifdef SPI   
  spi_HaloExchange<double> voltageExchange(sim.sendMap_, (sim.commTable_));
#else
  mpi_HaloExchange<double> voltageExchange(sim.sendMap_, (sim.commTable_));
#endif
  int nLocal = sim.anatomy_.nLocal();
  vector<double> dVmExternal(nLocal, 0.0);
  vector<double> dVmDiffusion(nLocal, 0.0);
  while ( sim.loop_<=sim.maxLoop_ )
  {
    int nLocal = sim.anatomy_.nLocal();
    
    //static TimerHandle haloHandle = profileGetHandle("Halo Exchange");
    //profileStart(haloHandle);
    voltageExchange.execute(sim.VmArray_, nLocal);
    voltageExchange.complete();
    //profileStop(haloHandle);

    for (unsigned ii=0; ii<nLocal; ++ii) dVmDiffusion[ii] = 0;
    for (unsigned ii=0; ii<nLocal; ++ii) dVmExternal[ii] = 0;
    

    // DIFFUSION
    //static TimerHandle diffusionHandle = profileGetHandle("Diffusion");
    //profileStart(diffusionHandle);
    sim.diffusion_->calc(sim.VmArray_, dVmDiffusion, voltageExchange.get_recv_buf_(), nLocal);
    //profileStop(diffusionHandle);

    //static TimerHandle stimulusHandle = profileGetHandle("Stimulus");
    //profileStart(stimulusHandle);
    for (unsigned ii=0; ii<sim.stimulus_.size(); ++ii) sim.stimulus_[ii]->stim(sim.time_, dVmDiffusion, dVmExternal);
    for (unsigned ii=0; ii<nLocal; ++ii) dVmDiffusion[ii] += dVmExternal[ii];
    //profileStop(stimulusHandle);
    L2_BarrierWithSync_WaitAndReset(&reactionBarrier, &reactionHandle, nReact);
    //static TimerHandle integratorHandle = profileGetHandle("Integrator");
    //profileStart(integratorHandle);
    for (unsigned ii=0; ii<nLocal; ++ii)
    {
      double dVm = dVmReaction[ii] + dVmDiffusion[ii];
      sim.VmArray_[ii] += sim.dt_*dVm;
    }
    sim.time_ += sim.dt_;
    ++sim.loop_;
    //profileStop(integratorHandle);

    L2_BarrierWithSync_Arrive(&diffusionBarrier, &diffusionHandle, nDiffuse);
    L2_BarrierWithSync_Reset(&diffusionBarrier, &diffusionHandle, nDiffuse);
    // SENSORS
    //static TimerHandle sensorHandle = profileGetHandle("Sensors");
    //profileStart(sensorHandle);
    for (unsigned ii=0; ii<sim.sensor_.size(); ++ii)
    {
      if (sim.loop_ % sim.sensor_[ii]->evalRate() == 0)
         sim.sensor_[ii]->eval(sim.time_, sim.loop_, sim.VmArray_, dVmReaction, dVmDiffusion, dVmExternal);
    }
    //profileStop(sensorHandle);
    for (unsigned ii=0; ii<nLocal; ++ii) dVmDiffusion[ii] -= dVmExternal[ii];
    loopIO(sim, dVmReaction, dVmDiffusion, dVmExternal, myRank);
    }
}
void reactionLoop(Simulate& sim, vector<double>& dVmReaction,L2_BarrierHandle_t& reactionHandle, L2_BarrierHandle_t& diffusionHandle)
{
  int ompTid = omp_get_thread_num();
  printf("reactionLoop: ompTid=%d\n",ompTid); 
  while ( sim.loop_<=sim.maxLoop_ )
  {
    int nLocal = sim.anatomy_.nLocal();
      
    //static TimerHandle reactionHandle = profileGetHandle("Reaction");
    //profileStart(reactionHandle);
    sim.reaction_->updateNonGate(sim.dt_, sim.VmArray_, dVmReaction);
    sim.reaction_->updateGate(sim.dt_, sim.VmArray_);
   // profileStop(reactionHandle);
     L2_BarrierWithSync_Arrive(&reactionBarrier, &reactionHandle, nReact);
     L2_BarrierWithSync_Reset(&reactionBarrier, &reactionHandle, nReact);
     L2_BarrierWithSync_WaitAndReset(&diffusionBarrier, &diffusionHandle, nDiffuse);
  }

}

void simulationLoopParallelDiffusionReaction(Simulate& sim)
{
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  int nLocal = sim.anatomy_.nLocal();
  vector<double> dVmReaction(nLocal, 0.0);
  vector<double> dVmDiffusion(nLocal, 0.0);
  vector<double> dVmExternal(sim.anatomy_.nLocal(), 0.0);
  simulationProlog(sim, dVmReaction,dVmDiffusion,dVmExternal,myRank);
  #pragma omp parallel
  {
  int ompTid = omp_get_thread_num();

  #pragma omp master
  nReact = omp_get_num_threads() - nDiffuse;

  L2_BarrierHandle_t reactionHandle;
  L2_BarrierHandle_t diffusionHandle;
  L2_BarrierWithSync_InitInThread(&reactionBarrier, &reactionHandle);
  L2_BarrierWithSync_InitInThread(&diffusionBarrier, &diffusionHandle);

  #pragma omp barrier
  if ( ompTid < nDiffuse ) diffusionLoop(sim, dVmReaction,  myRank, reactionHandle, diffusionHandle);
  else reactionLoop(sim, dVmReaction, reactionHandle, diffusionHandle);
  }
}
