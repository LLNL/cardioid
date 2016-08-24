
#include <mpi.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>

#include "Reaction.hh"
#include "VectorDouble32.hh"
#include "reactionFactory.hh"
#include "Anatomy.hh"
#include "ThreadServer.hh"

#include "singleCellOptions.h"

MPI_Comm COMM_LOCAL = MPI_COMM_WORLD;

/*!
  Keeps track of real time versus simulation ticks.
 */
class Timeline {
 public:
  Timeline(double dt, double duration) {
    maxTimesteps_ = round(duration/dt);
    dt_ = duration/maxTimesteps_;
  }
  int maxTimesteps() const { return maxTimesteps_; };
  double dt() const { return dt_; }
  double maxTime() const { return dt_*maxTimesteps_; }
  double realTimeFromTimestep(int timestep) const {
    return timestep*dt_;
  }
  int timestepFromRealTime(double realTime) const {
    return round(realTime/dt_);
  }
  
 private:
  double dt_;
  int maxTimesteps_;
};


int main(int argc, char* argv[]) {

  int npes, mype;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);  

  struct gengetopt_args_info params;
  cmdline_parser(argc, argv, &params);

  double dt = params.dt_arg;
  
  double duration;
  if (params.duration_given) {
    duration = params.duration_arg;
  } else {
    duration = params.s1_offset_arg + params.s1_count_arg*params.s1_bcl_arg;
  }

  Timeline timeline(dt, duration);


  //Stimulus setup
  std::vector<int> stimTimesteps;
  if (params.stim_at_given) {
    for(int istim=0; istim<params.stim_at_max; istim++) {
      stimTimesteps.push_back(timeline.timestepFromRealTime(params.stim_at_arg[istim]));
    }
  }
  bool s1_set_on_command_line = params.s1_count_given || params.s1_bcl_given || params.s1_offset_given; 
  if (s1_set_on_command_line || !params.stim_at_given) {
    for (int ibeat=0; ibeat<params.s1_count_arg; ibeat++) {
      stimTimesteps.push_back(timeline.timestepFromRealTime
                              (params.s1_offset_arg + ibeat*params.s1_bcl_arg));
    }
  }
  std::sort(stimTimesteps.begin(), stimTimesteps.end());

  int timestepsInEachStimulus=timeline.timestepFromRealTime(params.stim_duration_arg);
  double stimulusStrength = params.stim_strength_arg;
  int timestepsLeftInStimulus=0;
  std::vector<int>::iterator nextStimulus=stimTimesteps.begin();


  //Output setup
  int outputTimestepInterval = timeline.timestepFromRealTime(params.output_dt_arg);


  
  //create the ionic model
  const int nCells = 1;
  std::string modelName = params.model_arg;
  Anatomy anatomy;
  {
    anatomy.setGridSize(nCells,1,1);
    std::vector<AnatomyCell>& cells = anatomy.cellArray();
    cells.reserve(nCells);
    for (int ii=0; ii<nCells; ii++) {
      AnatomyCell tmp;
      tmp.gid_ = ii;
      tmp.cellType_ = 102;
      cells.push_back(tmp);
    }
  }
  ThreadServer& threadServer = ThreadServer::getInstance();
  ThreadTeam threads = threadServer.getThreadTeam(std::vector<unsigned>());
  Reaction* reaction = reactionFactory(modelName,dt,
                                       anatomy, threads,
                                       std::vector<std::string>());

  //initialize the ionic model
  VectorDouble32 Vm(nCells);
  std::vector<double> iStim(nCells);
  VectorDouble32 dVm(nCells);
  reaction->initializeMembraneVoltage(Vm);

  //read from a checkpoint if necessary.


  //using a forever loop here so we can repeat outputs on the last iteration.
  //otherwise, we'd use a for loop here.
  int itime=0;
  while(1) {
    //if we should checkpoint, do so.

    if (itime % outputTimestepInterval == 0) {
      //doIO();
      printf("%.16g %.16g\n", timeline.realTimeFromTimestep(itime), Vm[0]);
    }
  
    //if we're done with the loop, exit.
    if (itime == timeline.maxTimesteps()) { break; }

    //any stimulii left?
    if (nextStimulus != stimTimesteps.end()) {
      //if so, wait for it.
      if (itime == *nextStimulus) {
        ++nextStimulus;
        timestepsLeftInStimulus=timestepsInEachStimulus;
      }
    }
    for (int ii=0; ii<nCells; ii++) {
      iStim[ii]=0;
      if (timestepsLeftInStimulus > 0) {
        /* Check the negative sign here.  Look at:

           startTimer(stimulusTimer);
           // add stimulus to dVmDiffusion
           for (unsigned ii=0; ii<sim.stimulus_.size(); ++ii)
              sim.stimulus_[ii]->stim(sim.time_, dVmDiffusion);
           for (unsigned ii=0; ii<nLocal; ++ii)
              iStim[ii] = -(dVmDiffusion[ii]);
           stopTimer(stimulusTimer);

           Notice that the stimulii all add their current to
           stimulus. If you assume diffusion is zero, then istim ends
           up negative.

        */
        iStim[ii] = -params.stim_strength_arg;
        timestepsLeftInStimulus--;
      }
    }

    //calc() dVm and update the internal states
    reaction->calc(dt, Vm, iStim, dVm);
    
    //update Vm and apply stimulus
    for (int ii=0; ii<nCells; ii++) {
      //use a negative sign here to undo the negative we had above.
      Vm[ii] += (dVm[ii] - iStim[ii]) * dt;
    }

    itime++;
  }


  MPI_Finalize();
  return 0;
}
