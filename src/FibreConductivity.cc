#include "FibreConductivity.hh"
#include "AnatomyCell.hh"
#include <cmath>
#include <iostream>

inline void
calcConductivityMatrixIBT(SymmetricTensor& conductivity,
                          double sigmal, double sigmat,
                          double phi, double theta);



FibreConductivity::FibreConductivity(const FibreConductivityParms& p)
: sigmaTi_(p.sigmaTi), sigmaLi_(p.sigmaLi)
{
}


void FibreConductivity::compute(double theta, double phi, SymmetricTensor& sigma)
{
   calcConductivityMatrixIBT(sigma, sigmaLi_, sigmaTi_, phi, theta);
}


/** Adapted from conductivity.h in BlueBeats code.
 */
inline void
calcConductivityMatrixIBT(SymmetricTensor& conductivity,
                          double sigmal, double sigmat,
                          double phi, double theta)
{
   double sinPhi   = std::sin(  phi); // SinTable[phi];
   double cosPhi   = std::cos(  phi); // CosTable[phi];
   double sinTheta = std::sin(theta); // SinTable[theta];
   double cosTheta = std::cos(theta); // CosTable[theta];
   
   
   double sinPhiSquared = sinPhi * sinPhi;
   double cosPhiSquared = cosPhi * cosPhi;
   double sinThetaSquared = sinTheta * sinTheta;
   double cosThetaSquared = cosTheta * cosTheta;
   
   // Equations can still be optimized computationally!!!
   
   //M(0,0)
   conductivity.a11 = ((sigmal) * cosPhiSquared * sinThetaSquared) +
      ((sigmat) * ((cosPhiSquared * cosThetaSquared) + sinPhiSquared));
   //M(1,1)
   conductivity.a22 = ((sigmal) * sinPhiSquared * sinThetaSquared)+
      ((sigmat) * (cosPhiSquared + (sinPhiSquared*cosThetaSquared)));
   //M(2,2)
   conductivity.a33 = ((sigmal) * cosThetaSquared)+
      ((sigmat) *  sinThetaSquared);
   //M(0,1) = M(1,0)
   conductivity.a12 = ((sigmal)  -
                       (sigmat)) * cosPhi * sinThetaSquared * sinPhi;
   //M(0,2) = M(2,0)
   conductivity.a13 = ((sigmal) -
                       (sigmat)) * cosTheta * cosPhi * sinTheta;
   //M(1,2) = M(2,1)
   conductivity.a23 = ((sigmal) -
                       (sigmat)) * sinPhi * cosTheta * sinTheta;

   //ewd
   if (false && conductivity.a11 <= 0.0)
   {
      std::cout
         << sinPhi << " " << cosPhi << " " << sinTheta << " " << cosTheta << "\n"
         << conductivity.a11 << "\t"
         << conductivity.a12 << "\t"
         << conductivity.a13 << "\n"
         << conductivity.a22 << "\t"
         << conductivity.a23 << "\n"
         << conductivity.a33 << "\t"
         << std::endl;
   }                         

}
