#ifndef CLOOPER_H
#define CLOOPER_H

#ifdef __cplusplus
extern "C" {
#endif

#include <inttypes.h>
/** dVmR and dVmD really should be passed in as const double* since they
 * aren't modified.  However, I think this might interfere with
 * simdization */
/* void integrateLoop(const int begin, const int end, const double dt, */
/*                    double* dVmR, double* dVmD, double* Vm); */

void integrateLoop(const int begin, const int end,
                   const double dt,
                   double* dVmR,
                   double* stim,
                   unsigned* blockOffset,
                   double* dVmBlock,
                   double* VmBlock,
                   double* Vm,
                   double diffusionScale);
void integrateLoop_nostim(const int begin, const int end,
                   const double dt,
                   double* dVmR,
                   double* stim,
                   unsigned* blockOffset,
                   double* dVmBlock,
                   double* VmBlock,
                   double* Vm,
                   double diffusionScale);
void integrateLoopLag(const int begin, const int end,
                   const double dt,
                   double* dVmR,
                   double* stim,
                   unsigned* blockOffset,
                   double* dVmBlock,
                   double* VmBlock,
                   double* Vm,
                   double diffusionScale);


#ifdef __cplusplus
} // extern "C"
#endif

#endif
