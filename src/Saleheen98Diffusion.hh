#ifndef SALEHEEN98_DIFFUSION_HH
#define SALEHEEN98_DIFFUSION_HH

#include "Diffusion.hh"
#include <string>
#include "Array3d.hh"
#include "LocalGrid.hh"
#include <inttypes.h>
#include "ibmIntrinsics.hh"

class Anatomy;
class SymmetricTensor;
class Tuple;



struct Saleheen98DiffusionParms
{
   double diffusionScale_;
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
class Saleheen98PrecomputeDiffusion : public Diffusion
{
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
   
 public:
   Saleheen98PrecomputeDiffusion(
      const Saleheen98DiffusionParms& parms,
      const Anatomy& anatomy);

   
   void calc(const std::vector<double>& Vm, std::vector<double>& dVm, double *recv_buf_, int nLocal);
   void calc_simd(const std::vector<double>& Vm, std::vector<double>& dVm, double *recv_buf_, int nLocal);
   
 private:
   void   buildTupleArray(const Anatomy& anatomy);
   void   buildBlockIndex(const Anatomy& anatomy);
   void   precomputeCoefficients(const Anatomy& anatomy);
   void   reorder_Coeff();  //reoder the diffusion coefficients for simdizaed calc
   
   double boundaryFDLaplacianSaleheen98SumPhi(const Tuple& tuple);
   void   boundaryFDLaplacianSaleheen98SumPhi_All_simd(const uint32_t start,const int32_t chunk_size, Array3d<double>*  VmTmp, double* out);
   void   boundaryFDLaplacianSaleheen98Constants(
      const int*** tissue,
      const SymmetricTensor*** sigmaMatrix,
      const int& x, const int& y, const int& z,
      const double& dxInv, const double& dyInv, const double& dzInv);
      
   
   void updateVoltageBlock(const std::vector<double>& Vm, double *recv_buf, int nLocal);
   void printAllConductivities(const Array3d<int>& tissue,
                               const Array3d<SymmetricTensor>& sigma);
   void printAllDiffusionWeights(const Array3d<int>& tissue);


   LocalGrid                      localGrid_;
   double                         diffusionScale_;
   std::vector<unsigned>          blockIndex_; // for local and remote cells
   std::vector<Tuple>             localTuple_; // only for local cells
   Array3d<DiffusionCoefficients> diffIntra_;
   Array3d<double>                VmBlock_;
   Array3d<double>                diffCoefT2_;
   
   
};


#endif
