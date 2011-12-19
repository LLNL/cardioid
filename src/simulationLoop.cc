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
    sim.tmap_["Halo Exchange"].start();
    voltageExchange.execute(sim.VmArray_, nLocal);
    sim.tmap_["Halo Exchange"].stop();

    for (unsigned ii=0; ii<nLocal; ++ii)
       dVmExternal[ii] = 0;
    
    // DIFFUSION
    sim.tmap_["Diffusion"].start();
    sim.diffusion_->calc(sim.VmArray_, dVmDiffusion);
    sim.tmap_["Diffusion"].stop();

    // code to limit or set iStimArray goes here.
    for (unsigned ii=0; ii<sim.stimulus_.size(); ++ii)
      sim.stimulus_[ii]->stim(sim.time_, dVmDiffusion, dVmExternal);

    for (unsigned ii=0; ii<nLocal; ++ii)
      iStim[ii] = -(dVmDiffusion[ii] + dVmExternal[ii]);
      
    // REACTION
    sim.tmap_["Reaction"].start();
    sim.reaction_->calc(sim.dt_, sim.VmArray_, iStim, dVmReaction);
    sim.tmap_["Reaction"].stop();
    for (unsigned ii=0; ii<nLocal; ++ii)
    {
      double dVm = dVmReaction[ii] + dVmDiffusion[ii] + dVmExternal[ii] ;
      sim.VmArray_[ii] += sim.dt_*dVm;
    }
    sim.time_ += sim.dt_;
    ++sim.loop_;

    // SENSORS
    for (unsigned ii=0; ii<sim.sensor_.size(); ++ii)
    {
      if (sim.loop_ % sim.sensor_[ii]->evalRate() == 0)
         sim.sensor_[ii]->eval(sim.time_, sim.loop_, sim.VmArray_,
                               dVmReaction, dVmDiffusion, dVmExternal);
      if (sim.loop_ % sim.sensor_[ii]->printRate() == 0)
         sim.sensor_[ii]->print(sim.time_, sim.loop_, sim.VmArray_,
                                dVmReaction, dVmDiffusion, dVmExternal);
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
      DirTestCreate(fullname.c_str());
      fullname += "/anatomy";
      writeCells(sim, fullname.c_str());
    }
      
  }
}
