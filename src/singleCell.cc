
#include <mpi.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include "singleCellOptions.h"

MPI_Comm COMM_LOCAL = MPI_COMM_WORLD;

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

  struct gengetopt_args_info params;
  cmdline_parser(argc, argv, &params);

  double duration;
  if (params.duration_given) {
    duration = params.duration_arg;
  } else {
    duration = params.s1_offset_arg + params.s1_count_arg*params.s1_bcl_arg;
  }

  Timeline timeline(params.dt_arg, duration);


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
  
  for (int itime=0; itime<timeline.maxTimesteps(); itime++) {

    //any stimulii left?
    if (nextStimulus != stimTimesteps.end()) {
      //if so, wait for it.
      if (itime == *nextStimulus) {
        ++nextStimulus;
        timestepsLeftInStimulus=timestepsInEachStimulus;
      }
    }
    double iStim=0;
    if (timestepsLeftInStimulus > 0) {
      iStim = params.stim_strength_arg;
      timestepsLeftInStimulus--;
    }

    //apply stimulus
    
    //calc();
    if (itime % outputTimestepInterval == 0) {
      //doIO();
    }
    //checkpoint();
    
  }
  
  return 0;
}
