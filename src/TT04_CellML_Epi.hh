#ifndef TT04_CELLML_EPI_HH
#define TT04_CELLML_EPI_HH

#include "TT04_CellML.hh"

class TT04_CellML_Epi : public TT04_CellML
{
 public:
   TT04_CellML_Epi();
   double calc(double Vm, double iStim, double states[17],
               double rates[17], double algebraic[67]);
   double defaultState(unsigned ii){return defaultState_[ii];}
 private:
   static double constants_[46];
   static double defaultState_[17];
};

#endif
