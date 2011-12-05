#ifndef TT06_CELLML_MID_HH
#define TT06_CELLML_MID_HH

#include "TT06_CellML.hh"

class TT06_CellML_Mid : public TT06_CellML
{
 public:
   TT06_CellML_Mid();
   double calc(double Vm, double iStim, double states[19],
               double rates[19], double algebraic[70]);
   double defaultState(unsigned ii){return defaultState_[ii];}
 private:
   static double constants_[53];
   static double defaultState_[19];
};

#endif
