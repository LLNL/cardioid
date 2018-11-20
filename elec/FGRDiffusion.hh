#ifndef FGRDIFFUSION_HH
#define FGRDIFFUSION_HH

#include "Diffusion.hh"
#include "LocalGrid.hh"
#include "Array3d.hh"
#include "FGRUtils.hh"

class Anatomy;
class Vector;
class SymmetricTensor;
class ThreadTeam;
class L2_Barrier_t;
class L2_BarrierHandle_t;

class FGRDiffusion : public Diffusion
{
 public:
   FGRDiffusion(
      const FGRUtils::FGRDiffusionParms& parms,
      const Anatomy& anatomy,
      const ThreadTeam& threadInfo,
      const ThreadTeam& reactionThreadInfo);
   
   void updateLocalVoltage(ro_mgarray_ptr<double> VmLocal);
   void updateRemoteVoltage(ro_mgarray_ptr<double> VmRemote);
   void calc(rw_mgarray_ptr<double> dVm);
   unsigned* blockIndex(){return &blockIndex_[0];}
   double* VmBlock() {return VmBlock_.cBlock();}
   double* dVmBlock(){return dVmBlock_.cBlock();}

 private:
   void FGRDiff_simd(const uint32_t start,const int32_t chunk_size, Array3d<double>* VmTmp, double* out);
   void FGRDiff_simd_thread(const uint32_t start,const int32_t chunk_size, Array3d<double>* VmTmp, double* out);

   void buildTupleArray(const Anatomy& anatomy);
   void buildBlockIndex(const Anatomy& anatomy);
   void precomputeCoefficients(const Anatomy& anatomy);
   void precomputeCoefficientsPositive(const Anatomy & anatomy);
   void reorder_Coeff();

   void mkTissueArray(const Array3d<int>& tissueBlk, int ib, int* tissue);
   Vector f1(int ib, int iFace, const Vector& h,
             const Array3d<SymmetricTensor>& sigmaBlk);
   void printAllWeights(const Array3d<int>& tissue);

   int                             nLocal_;
   int                             offset_[19];
   int                             offsetMap_[27];
   int                             faceNbrOffset_[6];
   LocalGrid                       localGrid_;
   const ThreadTeam&               threadInfo_;
   const ThreadTeam&               reactionThreadInfo_;
   L2_Barrier_t*                   fgrBarrier_;
   std::vector<L2_BarrierHandle_t> barrierHandle_;
   std::vector<int>                localCopyOffset_;
   std::vector<int>                remoteCopyOffset_;
   std::vector<int>                threadOffset_;
   std::vector<int>                threadOffsetSimd_;
   std::vector<unsigned>           blockIndex_; // for local and remote cells
   std::vector<Tuple>              localTuple_; // only for local cells
   Array3d<double>                 A0_;
   Array3d<FGRUtils::DiffWeight>   weight_;
   Array3d<double>                 VmBlock_;
   Array3d<WeightType>             diffCoefT2_;
   Array3d<double>                 dVmBlock_;
};

#endif
