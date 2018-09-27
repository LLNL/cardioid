#include "Simulate.hh"
#include "Sensor.hh"

#include <cmath>
#include <mpi.h>
#include "Diffusion.hh"
#include <stdio.h>

//using std::isnan;
using std::vector;


/** Offsets version.  Be sure to maintain in parallel with OMP version*/
void Simulate::checkRanges(int begin, int end,
                           ro_mgarray_ptr<double> _Vm,
                           ro_mgarray_ptr<double> _dVmReaction,
                           ro_mgarray_ptr<double> _dVmDiffusion)
{
   ro_array_ptr<double> Vm = _Vm.useOn(CPU);
   ro_array_ptr<double> dVmReaction = _dVmReaction.useOn(CPU);
   ro_array_ptr<double> dVmDiffusion = _dVmDiffusion.useOn(CPU);
   //const double vMax =   60.;
   //const double vMin = -110.;
   const double vMax = checkRange_.vMax;
   const double vMin = checkRange_.vMin;
   for (unsigned ii=begin; ii<end; ++ii)
   {
      if ( Vm[ii] > vMax || Vm[ii] < vMin )
         outOfRange(ii, Vm[ii], dVmReaction[ii], dVmDiffusion[ii]);
   }
}

/** Omp version.  Be sure to maintain in parallel with offsets version.
 *  Don't call from parallel loop */
void Simulate::checkRanges(ro_mgarray_ptr<double> _Vm,
                           ro_mgarray_ptr<double> _dVmReaction,
                           ro_mgarray_ptr<double> _dVmDiffusion)
{
   ro_array_ptr<double> Vm = _Vm.useOn(CPU);
   ro_array_ptr<double> dVmReaction = _dVmReaction.useOn(CPU);
   ro_array_ptr<double> dVmDiffusion = _dVmDiffusion.useOn(CPU);
   int nLocal = anatomy_.nLocal();
   //const double vMax =   60.;
   //const double vMin = -110.;
   const double vMax = checkRange_.vMax;
   const double vMin = checkRange_.vMin;
   #pragma omp parallel for
   for (int ii=0; ii<nLocal; ++ii)
   {
      if ( Vm[ii] > vMax || Vm[ii] < vMin )
         outOfRange(ii, Vm[ii], dVmReaction[ii], dVmDiffusion[ii]);
   }
}

void Simulate::outOfRange(unsigned index, const double Vm, const double dVmr, const double dVmd)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank); 

#if 0 //FIXME!!!!!!
   /** This is awful.  Some diffusion classes don't store the results in
    *  array form, but keep them in an internal matrix.  We have to go
    *  fetch them. */
   if (diffusion_->dVmBlock() != 0)
   {
      unsigned bi = diffusion_->blockIndex()[index];
      dVmd += diffusion_->dVmBlock()[bi] * diffusion_->diffusionScale();
   }
#endif
   
   printf("WARNING: Voltage out of range: rank %d, index %d gid %llu "
          "loop = %d, V = %e, dVmd = %e, dVmr = %e\n",
          myRank, index, anatomy_.gid(index), loop_, Vm, dVmd, dVmr);
   fflush(stdout);
}

// check if any IO may be needed at this step
bool Simulate::checkIO(int loop) const
{
   if ( loop<0 )loop=loop_;
   
   if (loop > 0 && checkpointRate_ > 0 && loop % checkpointRate_ == 0)return true;
   
   for (unsigned ii=0; ii<sensor_.size(); ++ii)
   {
      if ( sensor_[ii]->checkEvalAtStep(loop) )
         return true;
      if ( sensor_[ii]->checkPrintAtStep(loop) )
         return true;
   }
   
   return false;
}

void Simulate::bufferReactionData()
{
   for(std::vector<Sensor*>::const_iterator is = sensor_.begin(); 
                                            is!= sensor_.end();
                                          ++is)
   {
      (*is)->bufferReactionData(loop_);
   }
}

void Simulate::bufferReactionData(const int begin, const int end)
{
   for(std::vector<Sensor*>::const_iterator is = sensor_.begin(); 
                                            is!= sensor_.end();
                                          ++is)
   {
      (*is)->bufferReactionData(begin, end, loop_);
   }
}
