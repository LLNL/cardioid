#ifndef JHU_CONDUCTIVITY_HH
#define JHU_CONDUCTIVITY_HH

#include "IndexToTuple.hh"

class AnatomyCell;

// ToDo:
//
// 1.  There are probably more parameters that could be exported to the
// input deck.  These include the homogeneous fiber angle, and the
// non-homogeneous epiAngle and endoAngle.


/** This code is based on functions and compile time macros found in the
 *  BlueBeats code in conductivity.h and cardiacSetupMPI.h.
 *
 *  When BlueBeats.cpp calls GetDiffusionTensorJHU it passes only two
 *  unique conductivity values: param.sigmaLi and param.sigmaTi.  These
 *  two parameters are mapped to GL_, GT_, and GN_ in the same way as in
 *  BlueBeats (GL_ = sigmaLi, GT_ = GN_ = sigmaTi).
 *
 *  The function GetDiffusionTensorJHU has two different expressions for
 *  E1, E2, and E3 depending on the value of the macro ROTATION_MATRIX.
 *  We reproduce that behavior with a run time parameter.
 *
 *  In BlueBeats the getFiberAngleJHU function can return either a
 *  homogeneous or a non-homogeneous value depending on compile time
 *  macros.  We default to the non-homogeneous case, but allow a run
 *  time switch.
 *
 *  In BlueBeats the sheetAngle is a compile time parameter.  We choose
 *  a default of 45 degrees.
 */
struct JHUConductivityParms
{
   bool homogeneousFiber;
   int rotatationMatrix;
   int transmuralAxis;

   // size of Anatomy to initialize IndexToTuple functor.
   int nx, ny, nz;
   
   double sigmaTi;
   double sigmaLi;

   double sheetAngle;
};


class JHUConductivity
{
 public:
   JHUConductivity(const JHUConductivityParms& p);
   void compute(AnatomyCell& cell);
   
 private:

   double getFiberAngle(const Tuple& globalCoord);
   double getImbriAngle(const Tuple& globalCoord);
   double getSheetAngle(const Tuple& globalCoord);

   bool homogeneousFiber_;
   int rotatationMatrix_;
   
   double GL_;
   double GT_;
   double GN_;

   double sheetAngle_;
   double homogeneousFiberAngle_;
   double endoAngle_;
   double twistPerCell_;

   IndexToTuple i2t_;
};

#endif
