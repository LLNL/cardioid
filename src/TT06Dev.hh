#ifndef TT06DEV_HH
#define TT06DEV_HH
#define SQ(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#include "TT06DevFit.hh"

void initState(double *STATES,int cellType);
void initCnst();
double computeUpdates(double dt, double Vm, double* STATES, int cellType);

class TT06Dev 
{
 public:
   TT06Dev(int);
   double calc(double dt, double Vm, double iStim, double states[nStateVar]);
   double defaultState(unsigned ii){return defaultState_[ii];}
 private:
   static double defaultState_[nStateVar];
   double c9; 
   double states_[nStateVar];
   int cellType_;
};

#endif
