#ifndef CONDUCTIVITY_HH
#define CONDUCTIVITY_HH

class AnatomyCell;

struct SigmaTensorMatrix
{
   double a11, a12, a13;
   double      a22, a23; // symmetric matrix. Thus a21 = a12
   double           a33; // symmetric matrix. Thus a31 = a13 and a32 = a23
};


class Conductivity
{
 public:
   virtual ~Conductivity(){};
   virtual void compute(const AnatomyCell& cell, SigmaTensorMatrix& sigma) = 0;
   virtual SigmaTensorMatrix defaultValue() = 0;
};


#endif


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
 *  This code has a number of internal switches that are set at compile
 *  time.  I've attempted to recreate all of those switches with runtime
 *  parmaters, but there isn't quite enough information in the BlueBeats
 *  code to provide adequate documentation.
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

