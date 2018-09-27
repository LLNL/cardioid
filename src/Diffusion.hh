#ifndef DIFFUSION_HH
#define DIFFUSION_HH

#include "lazy_array.hh"

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
   inline Diffusion(double newDiffusionScale=0) : diffusionScale_(newDiffusionScale) {}
   virtual ~Diffusion(){};
   virtual void updateLocalVoltage(ro_mgarray_ptr<double> VmLocal) = 0;
   virtual void updateRemoteVoltage(ro_mgarray_ptr<double> VmRemote) = 0;
   virtual void calc(rw_mgarray_ptr<double> dVm) = 0;
   virtual void calc_overlap(rw_mgarray_ptr<double> dVm) {};
   virtual unsigned* blockIndex(){return 0;}
   virtual double* VmBlock(){return 0;}
   virtual double* dVmBlock(){return 0;}
   virtual double diffusionScale(){return diffusionScale_;}
   virtual void setDiffusionScale(double newDiffusionScale) { diffusionScale_ = newDiffusionScale; }
   virtual void  dump_VmBlock(int tmp){;}
   virtual void test() {return;};

 protected:
   double diffusionScale_;
   
};

#endif
