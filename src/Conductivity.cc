// ToDo:
//
// 1.  Inquire about the units for angles in CalcConductivityMatrixIBT.
// At the very least I think the divisor needs to be 254 instead of
// 256. 


/** The fibre orientation is computed by a long block of code starting
 *  at approximately line 1165 of BlueBeats.cpp.
 *
 *  There appear to be three possible cases:
 *  1.  Code reads theta and phi from data files.
 *  2.  Data is synthesized using a model. 
 *  3.  Hard coded values are set (#define TEST_CONDUCTIVITY)
 *
 *  In any case, the effect of this code is to populate sigmaMintra with
 *  values.
 *
 *  Memory for sigmatensorMatrix*** sigmaMintra is allocated at line 750
 *  Initial default values (theta=0, phi=0) are computed by calling
 *  CalcConductivityMatrixIBT at line 823.
 *
 *  The following external parameters are used by this calculation:
 *
 *  - sigmaLi
 *  - sigmaTi
 *
 *
 *  DATA READ FROM FILES
 *  ====================
 *
 *  The theta and phi data is read (possibly in chunks) and broadcast to
 *  all tasks.  Each task then loops over the entire chunk, looking for
 *  cells that are inside its local block and calling
 *  CalcConductivityMatrixIBT for such cells.
 *
 *  Special Notes:
 *
 *  1.  Theta and phi are read from disk as unsigned char.  Hence their
 *  values must be in the range [0, 255].
 *  2.  Theta == 255 or phi == 255 is treated as a special case.  In
 *  this case:
 *    sigmaMintra.a11 = 0.1;
 *    sigmaMintra.a22 = 0.1;
 *    sigmaMintra.a33 = 0.1;
 *    sigmaMintra.a12 = 0.0;
 *    sigmaMintra.a13 = 0.0;
 *    sigmaMintra.a23 = 0.0;
 *  3.  If any of the diagonal elements of sigmaMintra are negative or
 *  zero, special code is called to write the matrix.  This code appears
 *  to be buggy as is does not write the correct symmetric matrix
 *  elements in row 3.
 *
 *
 *
 *  DATA SYNTHESIZED CASE
 *  =====================
 *
 *  GetDiffusionTensorJHU is called for cell types 75-77, 100-102, and
 *  30-31.  Otherwise, sigmaMintra is unmodified from initial defaults.
 *
 *
 *
 *
 *  HARD CODED VALUES
 *  =================
 *  sigmaMintra.a11 = 0.024;
 *  sigmaMintra.a22 = 0.15;
 *  sigmaMintra.a33 = 0.024;
 *  sigmaMintra.a12 = 0.0;
 *  sigmaMintra.a13 = 0.0;
 *  sigmaMintra.a23 = 0.0;
 *
 *  These are only used when no data is read and the TEST_CONDUCTIVITY
 *  macro is set at compile time.



 *  

*/

/** conductivity.h contains:
 *
 *  structs: sigmatensor, sigmatensorMatrix, diffusion2D, diffusion
 *
 *  kaTabularizedFunction templated class (unused)
 *  TablularizeSin (unused)
 *  TablularizeCos (unused)
 *  CalcConductivityMatrixIBT
 *  getFiberAngleJHU (could be static)
 *  getImbriAngleJHU (could be static)
 *  getSheetAngleJHU (could be static)
 *  GetDiffusionTensorJHU 
 *
 *  Special Notes:  The macro SHEET_ANGLE used in getSheetAngleJHU is
 *  defined in cardiacSetupMPI.h
 */



/** Questions for IBM:
 *
 *  1.  What is the valid range for theta and phi?  From the code it
 *  appears to be [0, 254], However, the units are apparently not
 *  degrees.  Rather, in CalcConductivityMatrixIBT the angles are
 *  multiplied by PI/256.  This seems wrong for two reasons:  First the
 *  divisor should be 254 since the value 255 is specifically excluded
 *  in BlueBeats.cpp.  Second, since we multiply by PI we only have a
 *  range of angles from [0, 180].  This is ok for one of the angles,
 *  but does not allow the expression of all possible directions.  Is
 *  there a symmetry relation here that I'm missing?
 *  
*/

// This is according to the IBT transformation/rotation Matrix R as
// implemented in the code
/** Note that theta and phi are *not* in degrees.  */
inline int
calcConductivityMatrixIBT(SigmaTensorMatrix* conductivity,
			  double sigmal, double sigmat,
			  int phi, int theta)
{
   double sinPhi   = std::sin(  phi*PI/256); // SinTable[phi];
   double cosPhi   = std::cos(  phi*PI/256); // CosTable[phi];
   double sinTheta = std::sin(theta*PI/256); // SinTable[theta];
   double cosTheta = std::cos(theta*PI/256); // CosTable[theta];
                                                 
              
   double sinPhiSquared = sinPhi * sinPhi;
   double cosPhiSquared = cosPhi * cosPhi;
   double sinThetaSquared = sinTheta * sinTheta;
   double cosThetaSquared = cosTheta * cosTheta;
                      
   // Equations can still be optimized computationally!!!
                        
   //M(0,0)
   conductivity->a11 = ((sigmal) * cosPhiSquared * sinThetaSquared) +
      ((sigmat) * ((cosPhiSquared * cosThetaSquared) + sinPhiSquared));
   //M(1,1)
   conductivity->a22 = ((sigmal) * sinPhiSquared * sinThetaSquared)+
      ((sigmat) * (cosPhiSquared + (sinPhiSquared*cosThetaSquared)));
   //M(2,2)
   conductivity->a33 = ((sigmal) * cosThetaSquared)+
      ((sigmat) *  sinThetaSquared);
   //M(0,1) = M(1,0)
   conductivity->a12 = ((sigmal)  -
			(sigmat)) * cosPhi * sinThetaSquared * sinPhi;
   //M(0,2) = M(2,0)
   conductivity->a13 = ((sigmal) -
			(sigmat)) * cosTheta * cosPhi * sinTheta;
   //M(1,2) = M(2,1)
   conductivity->a23 = ((sigmal) -
			(sigmat)) * sinPhi * cosTheta * sinTheta;

   if (conductivity->a11 <= 0.0)
   {
      std::cout
	 << sinPhi << " " << cosPhi << " " << sinTheta << " " << cosTheta << "\n"
	 << conductivity->a11 << "\t"
	 << conductivity->a12 << "\t"
	 << conductivity->a13 << "\n"
	 << conductivity->a22 << "\t"
	 << conductivity->a23 << "\n"
	 << conductivity->a33 << "\t"
	 << std::endl;
   }                         

   return 0;
}
