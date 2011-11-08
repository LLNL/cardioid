#ifndef TT06_CELLML_HH
#define TT06_CELLML_HH

class TT06_CellML
{
 public:
   virtual ~TT06_CellML(){};
   virtual double calc(double dt, double Vm, double iStim) = 0;
};

#endif
