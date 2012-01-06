#ifndef TT06_RRG_HH
#define TT06_RRG_HH


class TT06_RRG
{
 public:
   TT06_RRG(int cellType);
   double calc(double dt, double Vm, double iStim);
   
 private:
   static double constants_[53];
   double states_[19];
};

#endif
