#include "TransportCoordinator.hh"

void setStimulus(OnDevice<ArrayView<double>> iStim,OnDevice<ConstArrayView<double>> dVmDiffusion);
void integrateVoltage(OnDevice<ArrayView<double>> vmarray,
                      OnDevice<ConstArrayView<double>> dVmReaction,
                      OnDevice<ConstArrayView<double>> dVmDiffusionRaw,
                      const double dt);
