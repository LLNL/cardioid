
#include "Ledger.hh"

__global__ void integrateVoltageKernel(double* Vm, const double* dVm, const double* iStim, const double dt, const int nCells)
{
   int ii = threadIdx.x + blockIdx.x*blockDim.x;
   if (ii >= nCells) { return; }
   //use a negative sign here to undo the negative we had above.
   Vm[ii] += (dVm[ii] - iStim[ii]) * dt;
}

void integrateVoltage(double* Vm, const double* dVm, const double* iStim, const double dt, const int nCells)
{
   int blockSize = 1024;
   
   integrateVoltageKernel<<<(nCells+blockSize-1)/blockSize,blockSize>>>(
      ledger_lookup(Vm),
      ledger_lookup(dVm),
      ledger_lookup(iStim),
      dt,
      nCells);
}
