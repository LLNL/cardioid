#ifndef REACTION_HH
#define REACTION_HH

#include <vector>

class Reaction
{
 public:
   virtual ~Reaction(){};
   virtual void calc(double dt,
                     const std::vector<double>& Vm,
                     const std::vector<double>& iStim,
                     std::vector<double>& dVm) = 0;
   /** Populates the Vm array with some sensible default initial
    * membrane voltage.  Vm will be the parallel to the local cells in
    * the anatomy that was used to create the concrete reaction class.
    */
   virtual void initializeMembraneVoltage(std::vector<double>& Vm) = 0;
};

#endif
