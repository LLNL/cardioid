#ifndef TT04_CELLML_EPI_HH
#define TT04_CELLML_EPI_HH

#include "TT04_CellML.hh"

class TT04_CellML_Epi : public TT04_CellML
{
 public:
   TT04_CellML_Epi();
   double calc(double dt, double Vm, double iStim);
 private:
   static double constants_[46];
   static double defaultState_[17];
   double states_[17];
};

#endif
