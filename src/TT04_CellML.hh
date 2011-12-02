#ifndef TT04_CELLML_HH
#define TT04_CELLML_HH

class TT04_CellML
{
 public:
   virtual ~TT04_CellML(){};
   virtual double calc(double Vm, double iStim, double states[17],
                       double rates[17], double algebraic[67]) = 0;
   virtual double defaultState(unsigned ii) = 0;
};

#endif
