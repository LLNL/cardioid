#ifndef TT06DEV_HH
#define TT06DEV_HH


class TT06Dev 
{
 public:
   TT06Dev(int);
   double calc(double dt, double Vm, double iStim, double states[19]);
   double defaultState(unsigned ii){return defaultState_[ii];}
 private:
   static double constants_[59];
   static double defaultState_[19];
   double states_[19];
   int cellType_;
};

#endif
