#include "ReactionFHN.hh"
#include "Anatomy.hh"
#include "object.h"
#include "reactionFactory.hh"

using namespace std;

ReactionFHN::ReactionFHN(const int numPoints)
: nCells_(numPoints)
{
   // Model Parameters
   double a = 0.13;
   b_ = 0.013;
   double c1 = 0.26;
   double c2 = 0.1;
   c3_ = 1.0;
   Vrest_ = -85.0;
   Vpeak_ = 40.0;
   double vamp = Vpeak_ - Vrest_; //125.0
   Vthresh_ = Vrest_ + a * vamp; //-68.75;

   double W = 1.0;
   W_.resize(nCells_, W);

   fhn1_ = c1/(vamp*vamp);
   fhn2_ = c2/vamp;
}



void ReactionFHN::calc(double dt, const VectorDouble32& Vm,
                       const vector<double>& iStim , VectorDouble32& dVm)
{
  for (unsigned ii=0; ii<nCells_; ++ii)
  {
     dVm[ii] = fhn1_ * ( (Vm[ii] - Vrest_) * (Vm[ii] - Vthresh_) * (Vpeak_ - Vm[ii]) ) - (fhn2_ * (Vm[ii] - Vrest_) * W_[ii]);
     double Vtmp = Vm[ii] + (dVm[ii] - iStim[ii]) * dt;
     W_[ii] += dt * b_ * (Vtmp - Vrest_ - c3_ * W_[ii]);
  }
}

void ReactionFHN::initializeMembraneVoltage(VectorDouble32& Vm)
{
   assert(Vm.size() >= nCells_);
   Vm.assign(Vm.size(), -85.0);
}

REACTION_FACTORY(FHN)(OBJECT*, const double, const int numPoints, const ThreadTeam&)
{
   // None of the FHN model parameters are currently wired to the
   // input deck.
   return new ReactionFHN(numPoints);
}
REACTION_FACTORY(FitzhughNagumo)(OBJECT* obj, const double dt, const int numPoints, const ThreadTeam& group)
{
   return reactionFactoryForFHN(obj,dt,numPoints,group);
}



/** Everything below this point is the code from IBM upon which the
 * implementation above is based. */
#if 0
//#include "stdio.h"


#include <IBM_FHN.h>


IBM_FHN::IBM_FHN()
{
  Init();
  
}; // IBM_FHN::IBM_FHN()
  
  
  
IBM_FHN::~IBM_FHN()
{
};
  
  
void IBM_FHN::Init()
{
  // Model Parameters
  a = 0.13;
  b = 0.013;
  c1 = 0.26;
  c2 = 0.1;
  c3 = 1.0;
  Cm = 0.01;
  vrest = -85.0;
  vpeak = 40.0;
  vamp = vpeak - vrest; //125.0
  vthresh = vrest + a * vamp; //-68.75;

  W = 1.0;                          
  iapp = 0.0;

};


double IBM_FHN::Calc(double tinc,  double VinV,  double i_external)
{
  dtms = tinc * 1000; // convert [s] to [ms]
  VmV = VinV * 1000; // convert [V] to [mV]
  iapp = i_external;
  this->ODE_RHS(dtms);
  return(VmV * 0.001); //switch back to V
          
}; //end of IBM_FHN::Calc()


//double FitzhughNagumoSundnesEtAl (double &V, double &W, double &dt, double &iapp)
int IBM_FHN::ODE_RHS(double dt)
{
  // This works for mV and ms - note that values need to be scaled otherwise

//  double Iapp = vamp * iapp;
//  Iapp = Iapp/Cm;           // Iapp = - (Istim/Cm) ???
    
// dy/dt = f(t, y)
// forward Euler method
// y(tn+1) = y(tn) + dt * f(tn, yn)

// Rogers & McCulloch:
// dv/dt = (c1 * v * (v - a) * (1 - v)) - (c2 * v * w) + iapp;
// with forward Euler method:
// v = v + (dt * ((c1 * v * (v - a) * (1 - v)) - (c2 * v * w) + iapp));


// dw/dt = b * (v - c3 * w);
// with forward Euler method:
// w = w + (dt * (b * (v - c3 * w)));

// Sundnes et al. introduced:
// vamp = vpeak - vrest; defined above
// V0 = vamp * v + vrest;
// W0 = vamp * w;


  dVdt = ((c1/(vamp*vamp))* (VmV - vrest) * (VmV - vthresh) * (vpeak - VmV)) - ((c2/vamp) * (VmV - vrest) * W) + iapp;
  VmV = VmV + dtms * (dVdt);
  // dW/dt = b * (V - vrest - c3 * W);
  W = W + dtms * (b * (VmV - vrest - (c3 * W)));


  // From Pullan
  // Equation 3.109
  //Iion = -Cm * (dVdt);
  //cout << " Iion " << Iion;
      
  // Eqution 3.131
  //Vm = Vm - (dt/Cm) * (Iion + Istim);
  //cout << " Vm " << Vm << endl;
        
  
  return (VmV);
};
#endif

