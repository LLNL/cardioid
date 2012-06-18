#include "Simulate.hh"
#include "Sensor.hh"

#include <cmath>
#include <mpi.h>
#include "Diffusion.hh"
#include <stdio.h>

using std::isnan;
using std::vector;


/** */
void Simulate::checkRanges(const vector<double>& dVmReaction,
                           const vector<double>& dVmDiffusion)
{
   unsigned nLocal = anatomy_.nLocal();
   for (unsigned ii=0; ii<nLocal; ++ii)
   {
      if ( vdata_.outOfRange(ii) )
         outOfRange(ii, (*vdata_.dVmReaction_)[ii], (*vdata_.dVmDiffusion_)[ii]);
   }
}

void Simulate::outOfRange(unsigned index, double dVmr, double dVmd)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank); 
   double Vm = (*vdata_.VmArray_)[index];

   /** This is awful.  Some diffusion classes don't store the results in
    *  array form, but keep them in an internal matrix.  We have to go
    *  fetch them. */
   if (diffusion_->VmBlock() != 0)
   {
      unsigned bi = diffusion_->blockIndex()[index];
      dVmd += diffusion_->VmBlock()[bi] * diffusion_->diffusionScale();
   }
   
   printf("WARNING: Voltage out of range: rank %d, index %d gid %llu\n"
          "         loop = %d, V = %e, dVmd = %e, dVmr = %e\n",
          myRank, index, anatomy_.gid(index), loop_, Vm, dVmd, dVmr);
   fflush(stdout);
}

// check if any IO may be needed at this step
bool Simulate::checkIO()const
{
   if( loop_ % printRate_ == 0 )return true;
   if( loop_ > 0 && checkpointRate_ > 0 
                 && loop_ % checkpointRate_ == 0)return true;
   if (loop_ > 0 && loop_ % snapshotRate_ == 0)return true;
   
   for(std::vector<Sensor*>::const_iterator is = sensor_.begin(); 
                                            is!= sensor_.end();
                                          ++is)
   {
      if( (*is)->checkPrintAtStep(loop_) )return true;
   }
   
   return false;
}
