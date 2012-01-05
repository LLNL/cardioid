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

using namespace std;

void simulationLoop(Simulate& sim)
{
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   
  vector<double> dVmDiffusion(sim.anatomy_.nLocal(), 0.0);
  vector<double> dVmReaction(sim.anatomy_.nLocal(), 0.0);
  vector<double> dVmExternal(sim.anatomy_.nLocal(), 0.0);
  vector<double> iStim(sim.anatomy_.nLocal(), 0.0);
   
  // for now, hardcode initialization of voltage.
  // use TT04 value from BlueBeats
  sim.VmArray_.resize(sim.anatomy_.size(), -85.23); // in mVolts
  for (unsigned ii=sim.anatomy_.nLocal(); ii<sim.anatomy_.size(); ++ii)
     sim.VmArray_[ii] = 0;
   
  HaloExchange<double> voltageExchange(sim.sendMap_, *(sim.commTable_));


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
    
    static TimerHandle barrierHandle = profileGetHandle("Barrier");
    profileStart(barrierHandle);
    MPI_Barrier(MPI_COMM_WORLD);
    profileStop(barrierHandle);

    static TimerHandle haloHandle = profileGetHandle("Halo Exchange");
    profileStart(haloHandle);
    voltageExchange.execute(sim.VmArray_, nLocal);
    profileStop(haloHandle);

    for (unsigned ii=0; ii<nLocal; ++ii)
       dVmExternal[ii] = 0;
    
    // DIFFUSION
    static TimerHandle diffusionHandle = profileGetHandle("Diffusion");
    profileStart(diffusionHandle);
    sim.diffusion_->calc(sim.VmArray_, dVmDiffusion);
    profileStop(diffusionHandle);

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
}
