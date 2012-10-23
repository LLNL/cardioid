#ifndef DIFFUSION_HH
#define DIFFUSION_HH

#include "VectorDouble32.hh"

/**
 *  We have managed to do something exceptionally stupid in this class.
 *  The behavior of the calc function depends on the context from which
 *  it is called.  From the omp simulation loop the calc function should
 *  set the value of dVm for all cells.  From the
 *  parallelDiffusionReaction loop calc should increment the value of
 *  dVm.
 *
 *  The very fact that we have some implementations that work on the omp
 *  loop and some for the other is another clear sign that something is
 *  wrong.  There probably should be some sort of property defined in
 *  the concrete classes to declare shich loop is supported.
 */
class Diffusion
{
 public:
   virtual ~Diffusion(){};
   virtual void updateLocalVoltage(const double* VmLocal) = 0;
   virtual void updateRemoteVoltage(const double* VmRemote) = 0;
   virtual void calc(VectorDouble32& dVm) = 0;
   virtual void calc_overlap(VectorDouble32& dVm) {};
   virtual unsigned* blockIndex(){return 0;}
   virtual double* VmBlock(){return 0;}
   virtual double* dVmBlock(){return 0;}
   virtual double diffusionScale(){return 1;}
   virtual void  dump_VmBlock(int tmp){;}
   virtual void test() {return;};
};

#endif
