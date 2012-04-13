#include "clooper.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <unistd.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>

// dfr: called for all builds.  Must be portable.
void integrateLoop(const int begin, const int end,
                   const double dt,
                   double* dVmR,
                   double* stim,
                   unsigned* blockOffset,
                   double* dVmBlock,
                   double* VmBlock,
                   double* Vm,
                   double diffusionScale)
{
   for (unsigned ii=begin; ii<end; ++ii)
   {
      int ib = blockOffset[ii];

      double dVm = dVmR[ii] + stim[ii] + diffusionScale * dVmBlock[ib];
      Vm[ii] += dt*dVm;
      VmBlock[ib] = Vm[ii];
   }
   return;
}

void integrateLoop2(const int size, const double dt, double* dVmR, double* dVmD, double* Vm)
{
   const int halfSize = size/2;
   for (unsigned ii=0; ii<halfSize; ++ii)
   {
      double dVm1 = dVmR[ii] + dVmD[ii];
      double dVm2 = dVmR[ii+halfSize] + dVmD[ii+halfSize];
      Vm[ii] += dt*dVm1;
      Vm[ii+halfSize] += dt*dVm2;
   }
   if (size%2 != 0)
   {
      const int jj = size-1;
      double dVm1 = dVmR[jj] + dVmD[jj];
      Vm[jj] += dt*dVm1;
   }
   return;
}
