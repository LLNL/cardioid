#include "clooper.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <unistd.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>

// ewd:  These functions are loops that the BG/Q C++ compiler can't SIMDize.
// ewd:  Currently, these are only called when -DBGQ is set, so feel free
// ewd:  to put machine-specific code in here.

void integrateLoop(const int size, const double dt, double* dVmR, double* dVmD, double* Vm)
{
   for (unsigned ii=0; ii<size; ++ii)
   {
      double dVm = dVmR[ii] + dVmD[ii];
      Vm[ii] += dt*dVm;
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
