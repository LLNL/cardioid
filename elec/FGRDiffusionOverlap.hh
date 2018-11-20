#ifndef FGRDIFFUSIONOVERLAP_HH
#define FGRDIFFUSIONOVERLAP_HH

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

class FGRDiffusionOverlap : public Diffusion
{
 public:
   FGRDiffusionOverlap(
      const FGRUtils::FGRDiffusionParms& parms,
      const Anatomy& anatomy,
      const ThreadTeam& threadInfo,
      const ThreadTeam& reactionThreadInfo);
   
   void updateLocalVoltage(ro_mgarray_ptr<double> VmLocal);
   void updateRemoteVoltage(ro_mgarray_ptr<double> VmRemote);
   void calc(rw_mgarray_ptr<double> dVm);
   void calc_overlap(rw_mgarray_ptr<double> dVm);
   void test();
   unsigned* blockIndex(){return &blockIndex_[0];}
   double* VmBlock() {return VmBlock_.cBlock();}
   double* dVmBlock(){return dVmBlock_.cBlock();}

 private:
   void FGRDiff_simd_thread(const uint32_t start,const int32_t chunk_size, Array3d<double>* VmTmp, double* out);
   void FGRDiff_strip(const uint32_t start,const int32_t chunk_size, Array3d<double>* VmTmp, double* out);
   void FGRDiff_2D_xy(uint32_t slabID,const uint32_t start,const int32_t ey);
   void FGRDiff_2D_z(uint32_t slabID,const uint32_t start,const int32_t ey, double* out);
   void FGRDiff_3D(const uint32_t start,const int32_t chunk_size, Array3d<double>* VmTmp, double* out);

   void buildTupleArray(const Anatomy& anatomy);
   void buildBlockIndex(const Anatomy& anatomy);
   void buildSlabIndex(const Anatomy& anatomy);
   void precomputeCoefficients(const Anatomy& anatomy);
   void precomputeCoefficientsPositive(const Anatomy & anatomy);
   void reorder_Coeff();
   void reset_Coeff();

   void mkTissueArray(const Array3d<int>& tissueBlk, int ib, int* tissue);
   Vector f1(int ib, int iFace, const Vector& h,
             const Array3d<SymmetricTensor>& sigmaBlk);
   void printAllWeights(const Array3d<int>& tissue);
   void printAllVoltage(Array3d<double>& Voltage,int map);
   void setAllVoltage(Array3d<double>& Voltage,int cut=0);
   void copySlabToBlock();
   void clearBlockBD();
   void compareVoltage(Array3d<double>& A,Array3d<double>& B);

   int                             nLocal_;
   int                             nRemote_;
   int                             offset_[19];
   int                             offsetMap_[27];
   int                             faceNbrOffset_[6];
   uint64_t                        BBzb;
   uint64_t                        BBze;
   LocalGrid                       localGrid_;
   const ThreadTeam&               threadInfo_;
   const ThreadTeam&               reactionThreadInfo_;
   L2_Barrier_t*                   fgrBarrier_;
   std::vector<L2_BarrierHandle_t> barrierHandle_;
   std::vector<int>                localCopyOffset_;
   std::vector<int>                remoteCopyOffset_;
   std::vector<int>                threadOffset_;
   std::vector<int>                threadOffsetSimd_;
   std::vector<int>                threadOffset2D_[3];

   std::vector<unsigned>           blockIndex_; // for local and remote cells
   std::vector<unsigned>           slabIndex_; // for local and remote cells
   std::vector<Tuple>              localTuple_; // only for local cells
   Array3d<double>                 A0_;
   Array3d<FGRUtils::DiffWeight>   weight_;
   Array3d<double>                 VmBlock_;
   Array3d<WeightType>             diffCoefT2_;
   Array3d<double>                 dVmBlock_;

   double*                         VmSlab_[6]; // 6 slabs
   double*                         VmSlabT_; 
   unsigned                        VmSlabY_[6];
   unsigned                        VmSlabZ_[6];
   Array3d<WeightType>             diffCoefT3_[6]; //3x3 stencil on slabs. WeightType is to be redefined
   Array3d<double>                 dVmSlab_[6];
   unsigned                        dVmBlock2Doffset_[6];
   unsigned                        dVmBlock2Djump_[6];

  /* Stuff for yz-strip version. */
  std::vector< std::vector<int> >  strip_begin,strip_end;

};

#endif
