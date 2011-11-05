#ifndef TT04_CELLML_HH
#define TT04_CELLML_HH

class TT04_CellML
{
 public:
   virtual ~TT04_CellML(){};
   virtual double calc(double dt, double Vm, double iStim) = 0;
};

#endif
