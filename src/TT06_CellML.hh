#ifndef TT06_CELLML_HH
#define TT06_CELLML_HH

class TT06_CellML
{
 public:
   virtual ~TT06_CellML(){};
   virtual double calc(double Vm, double iStim, double states[19],
                       double rates[19], double algebraic[70]) = 0;
   virtual double defaultState(unsigned ii) = 0;
};

#endif
