#include "GPUDiffusion.hh"


GPUDiffusion::GPUDiffusion(const Anatomy& anatomy, int simLoopType)
: simLoopType_(simLoopType),
  blockIndex_(anatomy.size(), 0)
{
}

void GPUDiffusion::updateLocalVoltage(const double* VmLocal)
{
}

void GPUDiffusion::updateRemoteVoltage(const double* VmRemote)
{
}

void GPUDiffusion::calc(VectorDouble32& dVm)
{
   if (simLoopType_ == 0)
      dVm.assign(dVm.size(), 0.0);
}

unsigned* GPUDiffusion::blockIndex() {return &blockIndex_[0];}
double* GPUDiffusion::VmBlock() {return NULL;}
double* GPUDiffusion::dVmBlock() {return NULL;}


