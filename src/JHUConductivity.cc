#include "JHUConductivity.hh"
#include <cmath>
#include <cassert>
#include "AnatomyCell.hh"

using namespace std;

JHUConductivity::JHUConductivity(const JHUConductivityParms& p)
: homogeneousFiber_(p.homogeneousFiber),
  rotatationMatrix_(p.rotatationMatrix),
  GL_(p.sigmaLi),
  GT_(p.sigmaTi),
  GN_(p.sigmaTi),
  sheetAngle_(p.sheetAngle),
  homogeneousFiberAngle_(-60),
  endoAngle_(60), // degrees
  i2t_(p.nx, p.ny, p.nz)
{
   double epiAngle = -60; //degrees
   double totalTwist = endoAngle_ - epiAngle;    // 120 degrees
   // YDIM is TM axis (in JHU code) defined by max y coord
   twistPerCell_ =  totalTwist/p.transmuralAxis;
}

/** Adapted from GetDiffusionTensorJHU in BlueBeats conductivity.h */
void JHUConductivity::compute(AnatomyCell& cell)
{
   //D = R* [G_L 0 0; 0 G_T 0; 0 0 G_N]*R' where R = [E1 E2 E3] = Euler Vecotr Rotation Matrix
   
   Tuple globalCoord = i2t_(cell.gid_);
   double FiberAngle = getFiberAngle(globalCoord);
   double ImbriAngle = getImbriAngle(globalCoord);
   double SheetAngle = getSheetAngle(globalCoord);
   
   // Convert from Degrees to radians
   double FibAng = FiberAngle*M_PI/180.;
   double ImbrAng = ImbriAngle*M_PI/180.;
   double ShetAng = SheetAngle*M_PI/180.;
   
   // !!!!!!!!!!!!!!!
   // NEEDS MORE WORK - comment by Tabish
   // !!!!!!!!!!!!!!!!
   
   
   double E1[3];
   double E3[3];
   double E2[3];
   if (rotatationMatrix_ == 0)
   {
      // Use Katherine's Rotation Matrix
      
      E1[0] = cos(FibAng);               // E_1x
      E1[1] = 0;                         // E_1y
      E1[2] = sin(FibAng);               // E_1z
      E2[0] = -sin(FibAng)*sin(ShetAng); // E_2x
      E2[1] = cos(ShetAng);              // E_2y
      E2[2] = cos(FibAng)*sin(ShetAng);  // E_2z
      E3[0] = sin(FibAng)*cos(ShetAng);  // E_3x
      E3[1] = sin(ShetAng);              // E_3y
      E3[2] = -cos(FibAng)*cos(ShetAng); // E_3z
      
   }
   else if (rotatationMatrix_ == 1)
   {
      E1[0] = cos(FibAng);               // E_1x
      E1[1] = 0;                         // E_1y
      E1[2] = -sin(FibAng);              // E_1z
      E2[0] = sin(FibAng)*sin(ShetAng);  // E_2x
      E2[1] = cos(ShetAng);              // E_2y
      E2[2] = cos(FibAng)*sin(ShetAng);  // E_2z
      E3[0] = sin(FibAng)*cos(ShetAng);  // E_3x
      E3[1] = -sin(ShetAng);             // E_3y
      E3[2] = cos(FibAng)*cos(ShetAng);  // E_3z
   }
   else
      assert(false);
   
   // Check: Determinant of Rotation matrix 'E' should be 1
   double determinant_E =
      E1[0]*E2[1]*E3[2] + E2[0]*E3[1]*E1[2] + E3[0]*E1[1]*E2[2] -
      E1[0]*E2[2]*E3[1] - E2[0]*E1[1]*E3[2] - E3[0]*E2[1]*E1[2];

   assert( abs(determinant_E - 1) < 1e-10 );
   
   double D[3][3];
   D[0][0] = GL_*pow(E1[0],2) + GT_*pow(E2[0],2) + GN_*pow(E3[0],2);
   D[1][0] = GL_*E1[0]*E1[1] + GT_*E2[0]*E2[1] + GN_*E3[0]*E3[1];
   D[2][0] = GL_*E1[0]*E1[2] + GT_*E2[0]*E2[2] + GN_*E3[0]*E3[2];
   D[0][1] = D[1][0];
   D[1][1] = GL_*pow(E1[1],2) + GT_*pow(E2[1],2) + GN_*pow(E3[1],2);
   D[2][1] = GL_*E1[1]*E1[2] + GT_*E2[1]*E2[2] + GN_*E3[1]*E3[2];
   D[0][2] = D[2][0];
   D[1][2] = D[2][1];
   D[2][2] = GL_*pow(E1[2],2) + GT_*pow(E2[2],2) + GN_*pow(E3[2],2);

   for (unsigned ii=0; ii<3; ++ii)
      for (unsigned jj=0; jj<3; ++jj)
         if (abs(D[ii][jj]) < 1e-12 * GT_)
            D[ii][jj] = 0;
   
   
   // Check: Determinant of Diffusion Tensor 'D' should be 1
   double determinant_D =
      D[0][0]*D[1][1]*D[2][2] + D[0][1]*D[1][2]*D[2][0] + D[0][2]*D[1][0]*D[2][1] -
      D[0][0]*D[1][2]*D[2][1] - D[0][2]*D[1][1]*D[2][0] - D[0][1]*D[1][0]*D[2][2];
//   assert( abs(determinant_D - 1) < 1e-10 );
   
   // setting the conductivity matrix of the IBM code
   cell.sigma_.a11 = D[0][0];
   cell.sigma_.a12 = D[1][0];
   cell.sigma_.a13 = D[2][0];
   cell.sigma_.a22 = D[1][1];
   cell.sigma_.a23 = D[1][2];
   cell.sigma_.a33 = D[2][2];
}


double JHUConductivity::getFiberAngle(const Tuple& globalCoord)
{
   double fiberAngle;
   
   if (homogeneousFiber_)
   {
      // Homogeneous Wedge
      // Constant fiber angle, i.e. fibers are parallel to each other
      fiberAngle = homogeneousFiberAngle_;
   }
   else
   {
      // Transmural wedge
      // Fibers twist linearly from +60 degrees at Endo
      // to -60 degrees at Epi. twist angle = 120 degrees
      // fiber twist = twist angle/fiber length
      fiberAngle = endoAngle_ - (twistPerCell_ * (globalCoord.y()+1));
   }
   
  return fiberAngle;
}


double JHUConductivity::getImbriAngle(const Tuple&)
{
  double ImbrAng = 0;
  return ImbrAng;
}
      
double JHUConductivity::getSheetAngle(const Tuple&)
{
  return sheetAngle_;
}
