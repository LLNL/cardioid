#ifndef SALHEEN98_DIFFUSION_HH
#define SALHEEN98_DIFFUSION_HH

#include "Diffusion.hh"
#include <string>
#include "Array3d.hh"
#include "LocalGrid.hh"

class Anatomy;
class Conductivity;
class SigmaTensorMatrix;
class Tuple;



struct Salheen98DiffusionParms
{
   double diffusionScale_;
   std::string conductivityName_;
};


struct DiffusionCoefficients
{
   double sumA;
   double A1;
   double A2;
   double A3;
   double A4;
   double A5;
   double A6;
   double A7;
   double A8;
   double A9;
   double A10;
   double A11;
   double A12;
   double A13;
   double A14;
   double A15;
   double A16;
   double A17;
   double A18;
};

/** This Diffusion class takes advantage of the fact that the
 *  conductivities are constant in time so all of the diffusion
 *  coefficients can be precomputed.  The actual diffusion calculation
 *  is exactly what is done by the boundaryFDLaplacianSaleheen98SumPhi
 *  function call from BlueBeats.
 *
 *  This class takes care of all of the precomputation including:
 *
 *  - determines the size of the local grid (the bounding box for all
 *    local and remote cells in the Anatomy passed to the constructor)
 *  - builds and owns the array of Tuples (local grid indices) for all
 *    local and remote Anatomy cells.
 *  - Computes the conductivity matrix for all local and remote cells
 *  - Computes and owns the 19 constants in the DiffusionCoefficients
 *    structure for each local Anatomy cell.
 */
class Salheen98PrecomputeDiffusion : public Diffusion
{
 public:
   Salheen98PrecomputeDiffusion(
      const Salheen98DiffusionParms& parms,
      const Anatomy& anatomy);
   
   void calc(const std::vector<double>& Vm, std::vector<double>& dVm);
   
 private:
   void   buildTupleArray(const Anatomy& anatomy);
   void   buildBlockIndex(const Anatomy& anatomy);
   void   precomputeCoefficients(const Anatomy& anatomy);
   
   double boundaryFDLaplacianSaleheen98SumPhi(const Tuple& tuple);
   void   boundaryFDLaplacianSaleheen98Constants(
      const int*** tissue,
      const SigmaTensorMatrix*** sigmaMatrix,
      const int& x, const int& y, const int& z,
      const double& dxInv, const double& dyInv, const double& dzInv);
      
   
   void updateVoltageBlock(const std::vector<double>& Vm);


   LocalGrid                      localGrid_;
   double                         diffusionScale_;
   Conductivity*                  conductivity_;
   std::vector<unsigned>          blockIndex_; // for local and remote cells
   std::vector<Tuple>             localTuple_; // only for local cells
   Array3d<DiffusionCoefficients> diffIntra_;
   Array3d<double>                VmBlock_;
   
   
};


#endif
