#ifndef FIBRE_CONDUCTIVTY_HH
#define FIBRE_CONDUCTIVTY_HH

#include "SymmetricTensor.hh"

struct FibreConductivityParms
{
   double sigmaTi;
   double sigmaLi;
};


class FibreConductivity
{
 public:
   FibreConductivity(const FibreConductivityParms& p);
   
   void compute(double theta, double phi, SymmetricTensor& sigma);
   
 private:
   double sigmaTi_;
   double sigmaLi_;
};



#endif
