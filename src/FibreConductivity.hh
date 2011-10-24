#ifndef FIBRE_CONDUCTIVTY_HH
#define FIBRE_CONDUCTIVTY_HH

#include "Conductivity.hh"

struct FibreConductivityParms
{
   double sigmaTi;
   double sigmaLi;
};


class FibreConductivity : public Conductivity
{
 public:
   FibreConductivity(FibreConductivityParms& p);
   
   void compute(int theta, int phi, SigmaTensorMatrix& sigma);
 private:
   double sigmaTi_;
   double sigmaLi_;
};



#endif
