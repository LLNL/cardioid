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
   FibreConductivity(const FibreConductivityParms& p);
   
   void compute(const AnatomyCell& cell, SigmaTensorMatrix& sigma);
 private:
   double sigmaTi_;
   double sigmaLi_;
};



#endif
