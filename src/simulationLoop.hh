#ifndef SIMULATION_LOOP_HH
#define SIMULATION_LOOP_HH

class Simulate;

void simulationLoop(Simulate& sim);
void simulationLoopParallelDiffusionReaction(Simulate& sim);
void simulationLoopAllSkate(Simulate& sim);
void simulationLoopLag(Simulate& sim);

#endif
