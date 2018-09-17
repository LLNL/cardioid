#include "lazy_array.hh"

void setStimulus(wo_larray_ptr<double> iStim,ro_larray_ptr<double> dVmDiffusion);
void integrateVoltage(rw_larray_ptr<double> vmarray,
                      ro_larray_ptr<double> dVmReaction,
                      ro_larray_ptr<double> dVmDiffusionRaw,
                      const double dt);
