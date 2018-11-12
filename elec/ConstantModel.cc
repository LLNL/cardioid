#include "ConstantModel.hh"
#include <cmath>
#include <cassert>
#include <iostream>

using namespace std;

ConstantModel::ConstantModel(const std::vector<double>& eta,
                             const SymmetricTensor& sigma1,
                             const SymmetricTensor& sigma2,
                             const SymmetricTensor& sigma3,
                             const double alpha,
                             const double beta,
                             const double gamma):
                                sigma1_(sigma1),
                                sigma2_(sigma2),
                                sigma3_(sigma3),
                                alpha_(alpha),
                                beta_(beta),
                                gamma_(gamma)
{
   eta_[0]=eta[0];
   eta_[1]=eta[1];
   eta_[2]=eta[2];
   
   assert( fabs(eta_[0])>1.e-15 || fabs(eta_[1])>1.e-15 || fabs(eta_[2])>1.e-15 );
   assert( alpha_==alpha_ );
   assert( beta_==beta_ );
   assert( gamma_==gamma_ );
}

/** returns dVm/dt for the reaction part only. */
double ConstantModel::calc(const double x, const double y, const double z)
{
   assert( x>0. );
   assert( y>0. );
   assert( z>0. );
   assert( x<M_PI );
   assert( y<M_PI );
   assert( z<M_PI );
   double fx=-cos(x)*(
       eta_[1]*sin(z)*(sigma2_.a11*alpha_+sigma2_.a12*beta_+sigma2_.a13*gamma_)
      +eta_[2]*sin(y)*(sigma3_.a11*alpha_+sigma3_.a12*beta_+sigma3_.a13*gamma_));
   double fy=-cos(y)*(
       eta_[0]*sin(z)*(sigma1_.a12*alpha_+sigma1_.a22*beta_+sigma1_.a23*gamma_)
      +eta_[2]*sin(x)*(sigma3_.a12*alpha_+sigma3_.a22*beta_+sigma3_.a23*gamma_));
   double fz=-cos(z)*(
       eta_[0]*sin(y)*(sigma1_.a13*alpha_+sigma1_.a23*beta_+sigma1_.a33*gamma_)
      +eta_[1]*sin(x)*(sigma2_.a13*alpha_+sigma2_.a23*beta_+sigma2_.a33*gamma_));
      
   assert( fx==fx );
   assert( fy==fy );
   assert( fz==fz );
   
   return 1.e-3*(fx+fy+fz); // factor 1.e-3 because we used a conductivity in [mS]
}
