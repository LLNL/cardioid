#ifndef TT04DEV_HH
#define TT04DEV_HH

#include "TT04_CellML.hh"

class TT04Dev : public TT04_CellML
{
 public:
   TT04Dev(int);
   double calc(double dt, double Vm, double iStim);
 private:
   static double constants_[52];
   static double defaultState_[17];
   double states_[17];
   int cellType_; 
};

#endif
