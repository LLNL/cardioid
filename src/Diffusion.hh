#ifndef DIFFUSION_HH
#define DIFFUSION_HH

#include <vector>

class Diffusion
{
 public:
   virtual void
   calc(const std::vector<double>& Vm, std::vector<double>& dVm) = 0;
};

/** Notes:
 * 
 *  1. There are probably many things going on in the Salheen98Diffusion
 *     concrete class that should be moved to the base class.  For
 *     example setting up and handling of the local bounding box
 *     and multiplying by diffusionscale seem generic.  This needs to
 *     be worked out when we build the next derived class.
 */

#endif
