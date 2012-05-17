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

class FGRDiffusionFlex : public Diffusion
{
 public:
   FGRDiffusionFlex(
      const FGRUtils::FGRDiffusionParms& parms,
      const Anatomy& anatomy,
      const ThreadTeam& threadInfo,
      const ThreadTeam& reactionThreadInfo);
   
   void updateLocalVoltage(const double* VmLocal);
   void updateRemoteVoltage(const double* VmRemote);
   void calc(std::vector<double>& dVm);
   unsigned* blockIndex(){return &blockIndex_[0];}
   double* VmBlock() {return VmBlock_.cBlock();}
   double* dVmBlock(){return dVmBlock_.cBlock();}
   double diffusionScale(){return diffusionScale_;}

   //test functions
   void dump_anatomy( const Anatomy& anatomy,int slice);
   void dump_quad_anatomy( const Anatomy& anatomy,int slice);
   void dump_array3d(Array3d<double>& Box,int slice);
   void VmBlockSet();
   void test_calc(int slice);

 private:
   void FGRDiff_simd_thread(const uint32_t b_quad,const int32_t e_quad, double* out);

   void buildTupleArray(const Anatomy& anatomy);
   void buildBlockIndex(const Anatomy& anatomy);
   void precomputeCoefficients(const Anatomy& anatomy);
   void reorder_Coeff();
   void buildOffset(const Anatomy& anatomy);

   void mkTissueArray(const Array3d<int>& tissueBlk, int ib, int* tissue);
   Vector f1(int ib, int iFace, const Vector& h,
             const Array3d<SymmetricTensor>& sigmaBlk);
   void printAllWeights(const Array3d<int>& tissue);

   int                             nLocal_;
   int                             offset_[19];
   int                             faceNbrOffset_[6];
   LocalGrid                       localGrid_;
   double                          diffusionScale_;
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
   Array3d<FGRUtils::DiffWeight>   weight_;
   Array3d<double>                 VmBlock_;
   Array3d<WeightType>             diffCoefT2_;
   Array3d<double>                 dVmBlock_;

   int                             nx;
   int                             ny;
   int                             nz;
   int                             nTotal_;
   int                             nCalc_;
   int                             nQuad_;
   Array3d<int>                    QuadBlock_;
   std::vector<WeightType>         diffCoefT3_;
   std::vector<uint16_t>           dOffset_;
   std::vector<unsigned>           inIndex_; 
   std::vector<unsigned>           outIndex_; 

};

#endif
