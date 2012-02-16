#include "simulationLoop.hh"

#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>

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

using namespace std;

#ifdef TIMING
extern "C"{ long long timebase();}
#endif

void simulationLoop(Simulate& sim)
{
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

#ifdef TIMING
  const double cycles_to_usec = 1.0/850.0;  // BG/P = 850 MHz
  long long treact = 0.0;
#endif
  
  vector<double> dVmDiffusion(sim.anatomy_.nLocal(), 0.0);
  vector<double> dVmReaction(sim.anatomy_.nLocal(), 0.0);
  vector<double> dVmExternal(sim.anatomy_.nLocal(), 0.0);
  vector<double> iStim(sim.anatomy_.nLocal(), 0.0);
   
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
     BucketOfBits stateData =
        loadAndDistributeState(sim.stateFilename_, sim.anatomy_);
     assert(stateData.nRecords() == sim.anatomy_.nLocal());
     sim.reaction_->loadState(stateData);
     unsigned vmIndex = stateData.getIndex("Vm");
     if (vmIndex != stateData.nFields())
        for (unsigned ii=0; ii<stateData.nRecords(); ++ii)
           stateData.getRecord(ii).getValue(vmIndex, sim.VmArray_[ii]);
  }
  
  
#ifdef SPI   
  spi_HaloExchange<double> voltageExchange(sim.sendMap_, (sim.commTable_));
#else
  mpi_HaloExchange<double> voltageExchange(sim.sendMap_, (sim.commTable_));
#endif
  

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
     sim.sensor_[ii]->print(sim.time_, sim.loop_, sim.VmArray_,
                            dVmReaction, dVmDiffusion, dVmExternal);
  }

  while ( sim.loop_<=sim.maxLoop_ )
  {
    int nLocal = sim.anatomy_.nLocal();
    
#ifdef TIMING
    static TimerHandle barrierHandle = profileGetHandle("Barrier1");
    profileStart(barrierHandle);
    MPI_Barrier(MPI_COMM_WORLD);
    profileStop(barrierHandle);
#endif

#ifdef SPI
//    MPI_Barrier(MPI_COMM_WORLD);
//    voltageExchange.barrier();
#endif
    
    static TimerHandle haloHandle = profileGetHandle("Halo Exchange");
    profileStart(haloHandle);
#ifdef TIMING
    long long t1a = timebase();
    long long t1b = timebase();
#endif
    voltageExchange.execute(sim.VmArray_, nLocal);
    voltageExchange.complete();
#ifdef TIMING
    long long t2a = timebase();
    if (myRank == 0)
       cout << "TIMING: iteration " << sim.loop_ << ", haloExchange time = " << (t2a-t1b) << ", halo+timer = " << (t2a-t1a) << endl;
#endif
    profileStop(haloHandle);

    for (unsigned ii=0; ii<nLocal; ++ii)
       dVmExternal[ii] = 0;
    
#ifdef TIMING
    static TimerHandle barrierHandle2 = profileGetHandle("Barrier2");
    profileStart(barrierHandle2);
    MPI_Barrier(MPI_COMM_WORLD);
    profileStop(barrierHandle2);
#endif
    
    // DIFFUSION
    static TimerHandle diffusionHandle = profileGetHandle("Diffusion");
    profileStart(diffusionHandle);
//    sim.diffusion_->calc(sim.VmArray_, dVmDiffusion);
    sim.diffusion_->calc(sim.VmArray_, dVmDiffusion, voltageExchange.get_recv_buf_(), nLocal);
    profileStop(diffusionHandle);

#ifdef TIMING
    static TimerHandle barrierHandle3 = profileGetHandle("Barrier3");
    profileStart(barrierHandle3);
    MPI_Barrier(MPI_COMM_WORLD);
    profileStop(barrierHandle3);
#endif
    
    // code to limit or set iStimArray goes here.
    static TimerHandle stimulusHandle = profileGetHandle("Stimulus");
    profileStart(stimulusHandle);
    for (unsigned ii=0; ii<sim.stimulus_.size(); ++ii)
      sim.stimulus_[ii]->stim(sim.time_, dVmDiffusion, dVmExternal);

    for (unsigned ii=0; ii<nLocal; ++ii)
      iStim[ii] = -(dVmDiffusion[ii] + dVmExternal[ii]);
    profileStop(stimulusHandle);
      
    // REACTION
    static TimerHandle reactionHandle = profileGetHandle("Reaction");
    profileStart(reactionHandle);
#ifdef TIMING
    long long t3a = timebase();
#endif
    sim.reaction_->calc(sim.dt_, sim.VmArray_, iStim, dVmReaction);
#ifdef TIMING
    long long t4a = timebase();
    treact += (t4a-t3a);
    cout << "TIMING: myRank = " << myRank << ", iteration " << sim.loop_ << ", reaction time = " << (t4a-t3a) << " cycles, " << (t4a-t3a)*cycles_to_usec << " usec" << endl;
#endif
    profileStop(reactionHandle);

#ifdef TIMING
    static TimerHandle barrierHandle4 = profileGetHandle("Barrier4");
    profileStart(barrierHandle4);
    MPI_Barrier(MPI_COMM_WORLD);
    profileStop(barrierHandle4);
#endif
    
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
         sim.sensor_[ii]->eval(sim.time_, sim.loop_, sim.VmArray_,
                               dVmReaction, dVmDiffusion, dVmExternal);
      if (sim.loop_ % sim.sensor_[ii]->printRate() == 0)
         sim.sensor_[ii]->print(sim.time_, sim.loop_, sim.VmArray_,
                                dVmReaction, dVmDiffusion, dVmExternal);
    }
    profileStop(sensorHandle);
    
    
    static TimerHandle loopIOHandle = profileGetHandle("LoopIO");
    profileStart(loopIOHandle);
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

#ifdef TIMING  
  cout << "TIMING: myRank = " << myRank << ", total reaction time = " << (double)treact << " cycles, " << setprecision(14) << (double)treact*cycles_to_usec << " usec" << endl;
#endif


}
