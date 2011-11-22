#include "FibreConductivity.hh"
#include "AnatomyCell.hh"
#include <cmath>
#include <iostream>

inline void
calcConductivityMatrixIBT(SigmaTensorMatrix& conductivity,
			  double sigmal, double sigmat,
			  int phi, int theta);



FibreConductivity::FibreConductivity(const FibreConductivityParms& p)
: sigmaTi_(p.sigmaTi), sigmaLi_(p.sigmaLi)
{
}


void FibreConductivity::compute(const AnatomyCell& cell, SigmaTensorMatrix& sigma)
{
   int theta = cell.theta_;
   int phi = cell.phi_;
   calcConductivityMatrixIBT(sigma, sigmaLi_, sigmaTi_, phi, theta);
}



/** Adapted from conductivity.h in BlueBeats code.
 */
inline void
calcConductivityMatrixIBT(SigmaTensorMatrix& conductivity,
			  double sigmal, double sigmat,
			  int phi, int theta)
{
   double sinPhi   = std::sin(  phi*M_PI/256); // SinTable[phi];
   double cosPhi   = std::cos(  phi*M_PI/256); // CosTable[phi];
   double sinTheta = std::sin(theta*M_PI/256); // SinTable[theta];
   double cosTheta = std::cos(theta*M_PI/256); // CosTable[theta];
   
   
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
   
   if (conductivity.a11 <= 0.0)
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
