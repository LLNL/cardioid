#ifndef TT06_CELLML_ENDO_HH
#define TT06_CELLML_ENDO_HH

#include "TT06_CellML.hh"

class TT06_CellML_Endo : public TT06_CellML
{
 public:
   TT06_CellML_Endo();
   double calc(double Vm, double iStim, double states[19],
               double rates[19], double algebraic[70]);
   double defaultState(unsigned ii){return defaultState_[ii];}
 private:
   static double constants_[53];
   static double defaultState_[19];
};

#endif
