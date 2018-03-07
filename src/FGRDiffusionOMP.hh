#ifndef FGRDIFFUSION_OMP_HH
#define FGRDIFFUSION_OMP_HH

#include "Diffusion.hh"
#include "LocalGrid.hh"
#include "Array3d.hh"
#include "FGRUtils.hh"

class Anatomy;
class Vector;
class SymmetricTensor;

class FGRDiffusionOMP : public Diffusion
{
 public:
   FGRDiffusionOMP(
      const FGRUtils::FGRDiffusionParms& parms,
      const Anatomy& anatomy);
   
   void updateLocalVoltage(const Managed<ArrayView<double>> VmLocal);
   void updateRemoteVoltage(const Managed<ArrayView<double>> VmRemote);
   void calc(Managed<ArrayView<double>> dVm);

 private:
   void buildTupleArray(const Anatomy& anatomy);
   void buildBlockIndex(const Anatomy& anatomy);
   void precomputeCoefficients(const Anatomy& anatomy);

   void mkTissueArray(const Array3d<int>& tissueBlk, int ib, int* tissue);
   Vector f1(int ib, int iFace, const Vector& h,
             const Array3d<SymmetricTensor>& sigmaBlk);
   void printAllWeights(const Array3d<int>& tissue);

   int                             nLocal_;
   int                             nRemote_;
   int                             offset_[19];
   int                             faceNbrOffset_[6];
   LocalGrid                       localGrid_;
   double                          diffusionScale_;
   std::vector<unsigned>           blockIndex_; // for local and remote cells
   std::vector<Tuple>              localTuple_; // only for local cells
   Array3d<double>                 A0_;
   Array3d<FGRUtils::DiffWeight>   weight_;
   Array3d<double>                 VmBlock_;
};

#endif
