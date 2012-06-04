#include "FGRDiffusionFlex.hh"
#include <cassert>
#include "DiffusionUtils.hh"
#include "Anatomy.hh"
#include "Vector.hh"
#include "fastBarrier.hh"
#include <algorithm>
#include <cstdio>
#include "ibmIntrinsics.hh"
#include "ThreadServer.hh"
#include "ThreadUtils.hh"
#include "PerformanceTimers.hh"
using namespace PerformanceTimers;

using namespace std;
using namespace FGRUtils;

#define U16_MAX 65000 //max number of things are 65000/8

FGRDiffusionFlex::FGRDiffusionFlex(const FGRDiffusionParms& parms,
                           const Anatomy& anatomy,
                           const ThreadTeam& threadInfo,
                           const ThreadTeam& reactionThreadInfo)
: nLocal_(anatomy.nLocal()), nTotal_(anatomy.size()),
  localGrid_(DiffusionUtils::findBoundingBox_simd(anatomy, parms.printBBox_)),
  threadInfo_(threadInfo),
  reactionThreadInfo_(reactionThreadInfo),
  diffusionScale_(parms.diffusionScale_)
{

   nx = localGrid_.nx();
   ny = localGrid_.ny();
   nz = localGrid_.nz();

   assert(nz%4 == 0); //must be aligned

   // This is a test
   for (unsigned ii=0; ii<anatomy.size(); ++ii)
   {
      Tuple globalTuple = anatomy.globalTuple(ii);
      Tuple ll = localGrid_.localTuple(globalTuple);
      assert(ll.x() >= 0 && ll.y() >= 0 && ll.z() >= 0);
      assert(ll.x() < nx && ll.y() < ny && ll.z() < nz);
   }
   // This has been a test

   mkOffsets(threadOffset_,     anatomy.nLocal(),  threadInfo_);
   mkOffsets(localCopyOffset_,  anatomy.nLocal(),  reactionThreadInfo_);
   mkOffsets(remoteCopyOffset_, anatomy.nRemote(), threadInfo_);

   barrierHandle_.resize(threadInfo.nThreads());
   fgrBarrier_ = L2_BarrierWithSync_InitShared();
   #pragma omp parallel
   {
      int tid = threadInfo.teamRank();
      if (tid >= 0)
         L2_BarrierWithSync_InitInThread(fgrBarrier_, &barrierHandle_[tid]);
   }         
   
   std::cout << "nx,ny,nz=" << nx <<" "<< ny <<" "<< nz << std::endl;

   //Dual struecture co-exist.
   //Array3d for book-keeping  : blockIndex, localTuple 
   //normal Vector for calculation : inIndex, outIndex

   weight_.resize(nx, ny, nz);
   VmBlock_.resize(nx, ny, nz);
   dVmBlock_.resize(nx,ny,nz);

   buildTupleArray(anatomy);
   buildBlockIndex(anatomy);

   int base = VmBlock_.tupleToIndex(1, 1, 1);
   offset_[ZZZ] = 0;
   offset_[PZZ] = VmBlock_.tupleToIndex(2, 1, 1) - base;
   offset_[ZPZ] = VmBlock_.tupleToIndex(1, 2, 1) - base;
   offset_[MZZ] = VmBlock_.tupleToIndex(0, 1, 1) - base;
   offset_[ZMZ] = VmBlock_.tupleToIndex(1, 0, 1) - base;
   offset_[PPZ] = VmBlock_.tupleToIndex(2, 2, 1) - base;
   offset_[MPZ] = VmBlock_.tupleToIndex(0, 2, 1) - base;
   offset_[MMZ] = VmBlock_.tupleToIndex(0, 0, 1) - base;
   offset_[PMZ] = VmBlock_.tupleToIndex(2, 0, 1) - base;
   offset_[ZZP] = VmBlock_.tupleToIndex(1, 1, 2) - base;
   offset_[ZZM] = VmBlock_.tupleToIndex(1, 1, 0) - base;
   offset_[ZPP] = VmBlock_.tupleToIndex(1, 2, 2) - base;
   offset_[ZPM] = VmBlock_.tupleToIndex(1, 2, 0) - base;
   offset_[ZMM] = VmBlock_.tupleToIndex(1, 0, 0) - base;
   offset_[ZMP] = VmBlock_.tupleToIndex(1, 0, 2) - base;
   offset_[PZP] = VmBlock_.tupleToIndex(2, 1, 2) - base;
   offset_[MZP] = VmBlock_.tupleToIndex(0, 1, 2) - base;
   offset_[MZM] = VmBlock_.tupleToIndex(0, 1, 0) - base;
   offset_[PZM] = VmBlock_.tupleToIndex(2, 1, 0) - base;

   faceNbrOffset_[0] = offset_[MZZ];
   faceNbrOffset_[1] = offset_[ZMZ];
   faceNbrOffset_[2] = offset_[PZZ];
   faceNbrOffset_[3] = offset_[ZPZ];
   faceNbrOffset_[4] = offset_[ZZM];
   faceNbrOffset_[5] = offset_[ZZP];


//   offsetsTest();
   precomputeCoefficients(anatomy);

   //this is test
   boundary_test(anatomy);

   reorder_Coeff();
   buildOffset(anatomy);

   cout << "buildOffset done. " << endl;

   //this is test
   {
     Array3d<double> tmp_Vm;
     tmp_Vm.resize(nx,ny,nz,0.0);
     srand(1234);
     
     //moving into Array3d
     for(int ii=0;ii<nTotal_;ii++)
     {
       tmp_Vm(blockIndex_[ii]) = VmBlock_(inIndex_[ii]) = rand();
     }

     basic_diffusion(weight_,tmp_Vm,dVmBlock_);

     double* out = new double[nCalc_*4];

     
     stencil(0,1,out);
     stencil(1,nCalc_,out);

     //comparison
     for(int ii=0;ii<nLocal_;ii++)
     {
       double a = dVmBlock_(blockIndex_[ii]);
       double b = out[outIndex_[ii]];
       double big = a>b ? a:b;
       double small = a>b ? b:a;
       double difference = big-small;
       if ( big !=0 )  difference /= fabs(big);
       if ( difference > 0.001)
       {
          std::cout << "no match :" << ii << std::endl;
          Tuple localTuple = localTuple_[ii];
          int x = localTuple.x();
          int y = localTuple.y();
          int z = localTuple.z();
          std::cout << "x,y,z=" << x << " " << y  << " " << z << std::endl;
          assert(false);
       }
     }

     delete [] out;
   }

   cout << "nQuad,nCalc=" << nCalc_ << "," << nQuad_ << endl;

   //simd thread offsets
   int chunkSize = nCalc_ / threadInfo.nThreads();
   int leftOver = nCalc_ % threadInfo.nThreads();
   threadOffsetSimd_.resize(threadInfo.nThreads()+1);
   threadOffsetSimd_[0]=0;
   for (int ii=0; ii<threadInfo.nThreads(); ++ii)
   {
      threadOffsetSimd_[ii+1] = threadOffsetSimd_[ii] + chunkSize;
      if (ii < leftOver)
         ++threadOffsetSimd_[ii+1];
   }
   assert(nCalc_ == threadOffsetSimd_[threadInfo.nThreads()] );

}

void FGRDiffusionFlex::updateLocalVoltage(const double* VmLocal)
{
   startTimer(FGR_ArrayLocal2MatrixTimer);
   int tid = reactionThreadInfo_.teamRank();
   unsigned begin = localCopyOffset_[tid];
   unsigned end   = localCopyOffset_[tid+1];
   for (unsigned ii=begin; ii<end; ++ii)
   {
      int index = inIndex_[ii];
      VmBlock_(index) = VmLocal[ii];
   }
   stopTimer(FGR_ArrayLocal2MatrixTimer);
}

void FGRDiffusionFlex::updateRemoteVoltage(const double* VmRemote)
{
   startTimer(FGR_ArrayRemote2MatrixTimer);
   int tid = threadInfo_.teamRank();
   unsigned begin = remoteCopyOffset_[tid];
   unsigned end   = remoteCopyOffset_[tid+1];
   unsigned* bb = &inIndex_[nLocal_];
   for (unsigned ii=begin; ii<end; ++ii)
   {
      int index = bb[ii];
      VmBlock_(index) = VmRemote[ii];
   }
   stopTimer(FGR_ArrayRemote2MatrixTimer);
}

/** threaded simd version */
void FGRDiffusionFlex::calc(vector<double>& dVm)
{
   int tid = threadInfo_.teamRank();

   startTimer(FGR_AlignCopyTimer);
   
   stopTimer(FGR_AlignCopyTimer);
   
   startTimer(FGR_StencilTimer);

   //uint32_t begin = VmTmp->tupleToIndex(threadOffsetSimd_[tid],1,0);
   //uint32_t end = VmTmp->tupleToIndex(threadOffsetSimd_[tid+1]-1,VmTmp->ny()-2,VmTmp->nz());
   //   printf("simd version:%d-%d\n",begin,end);
   //if (threadOffsetSimd_[tid] < threadOffsetSimd_[tid+1] )
   //   FGRDiff_simd_thread(begin,end-begin,VmTmp,dVmBlock_.cBlock());
   if (threadOffsetSimd_[tid] < threadOffsetSimd_[tid+1] )
      stencil(threadOffsetSimd_[tid] ,threadOffsetSimd_[tid+1],dVmBlock_.cBlock());
//      FGRDiff_simd_thread(threadOffsetSimd_[tid] ,threadOffsetSimd_[tid+1],dVmBlock_.cBlock());

   stopTimer(FGR_StencilTimer);

   startTimer(FGR_Barrier2Timer);
   L2_BarrierWithSync_Barrier(fgrBarrier_, &barrierHandle_[tid], threadInfo_.nThreads());
   stopTimer(FGR_Barrier2Timer);

   startTimer(FGR_Matrix2ArrayTimer);

//    if(VmBlock_.nz()%4 != 0) delete VmTmp;

//     int begin = threadOffset_[tid];
//     int end   = threadOffset_[tid+1];
//     double* dVmBlock_ptr = dVmBlock_.cBlock();
//     for (int ii=begin; ii<end; ++ii)
//     {
//       //dVm[ii] = dVmBlock_(localTuple_[ii].x(),localTuple_[ii].y(),localTuple_[ii].z());
//       dVm[ii] = dVmBlock_ptr[blockIndex_[ii]];
//       dVm[ii] *= diffusionScale_;
//     }
   stopTimer(FGR_Matrix2ArrayTimer);
}

/** We're building the localTuple array only for local cells.  We can't
 * do stencil operations on remote particles so we shouldn't need
 * tuples.  We can use block indices instead.
 */
void FGRDiffusionFlex::buildTupleArray(const Anatomy& anatomy)
{
   localTuple_.resize(anatomy.size(), Tuple(0,0,0));
   for (unsigned ii=0; ii<anatomy.size(); ++ii)
   {
      Tuple globalTuple = anatomy.globalTuple(ii);
      localTuple_[ii] = localGrid_.localTuple(globalTuple);
      assert(localTuple_[ii].x() > -1);
      assert(localTuple_[ii].y() > -1);
      assert(localTuple_[ii].z() > -1);
      assert(localTuple_[ii].x() < localGrid_.nx());
      assert(localTuple_[ii].y() < localGrid_.ny());
      assert(localTuple_[ii].z() < localGrid_.nz());

//      cout << "localTuple " << ii << ":" << localTuple_[ii].x() << " " << localTuple_[ii].y() << " " << localTuple_[ii].z() << endl;
   }
}

void FGRDiffusionFlex::buildBlockIndex(const Anatomy& anatomy)
{
   blockIndex_.resize(anatomy.size());
   for (unsigned ii=0; ii<anatomy.size(); ++ii)
   {
      Tuple globalTuple = anatomy.globalTuple(ii);
      Tuple ll = localGrid_.localTuple(globalTuple);
      blockIndex_[ii] = VmBlock_.tupleToIndex(ll.x(), ll.y(), ll.z());
   }
}

void FGRDiffusionFlex::precomputeCoefficients(const Anatomy& anatomy)
{
   unsigned nxGlobal = anatomy.nx();
   unsigned nyGlobal = anatomy.ny();
   unsigned nzGlobal = anatomy.nz();
   Vector hInv(1.0/anatomy.dx(), 1.0/anatomy.dy(), 1.0/anatomy.dz());
   Vector h(anatomy.dx(), anatomy.dy(), anatomy.dz());
   double gridCellVolume = h[0]*h[1]*h[2];
   
   SymmetricTensor sigmaZero = {0};
   Array3d<SymmetricTensor> sigmaBlk(nx, ny, nz, sigmaZero);
   Array3d<int> tissueBlk(nx, ny, nz, 0);

   const vector<AnatomyCell>& cell = anatomy.cellArray();
   for (unsigned ii=0; ii<anatomy.size(); ++ii)
   {
      unsigned ib = blockIndex_[ii];
      sigmaBlk(ib) = anatomy.conductivity(ii);
      tissueBlk(ib) = isTissue(anatomy.cellType(ii));
   }

   for (unsigned ii=0; ii<weight_.size(); ++ii)
      for (unsigned jj=0; jj<19; ++jj)
         weight_(ii).A[jj] = 0.0;

   for (unsigned iCell=0; iCell<anatomy.nLocal(); ++iCell)
   {
      unsigned ib = blockIndex_[iCell];
      int tissue[19] = {0};
      mkTissueArray(tissueBlk, ib, tissue);
      
      for (unsigned iFace=0; iFace<6; ++iFace)
      {
         unsigned faceNbrIndex = ib+faceNbrOffset_[iFace];
         if (tissueBlk(faceNbrIndex) == 0)
            continue;

         Vector sigmaTimesS = f1(ib, iFace, h, sigmaBlk)/gridCellVolume;
         double gradPhi[3][19] = {0};
         f2(iFace, tissue, gradPhi);
         
         for (unsigned ii=0; ii<19; ++ii)
            for (unsigned jj=0; jj<3; ++jj)
               weight_(ib).A[ii] += sigmaTimesS[jj] * gradPhi[jj][ii] * hInv[jj];
      }
      double sum = 0;
      for (unsigned ii=0; ii<19; ++ii)
         sum += weight_(ib).A[ii];
      assert(abs(sum) < weightSumTolerance);
   }
//   printAllWeights(tissueBlk);
}

void FGRDiffusionFlex::mkTissueArray(
   const Array3d<int>& tissueBlk, int ib, int* tissue)
{
   for (unsigned ii=0; ii<19; ++ii)
      tissue[ii] = tissueBlk(ib + offset_[ii]);
}


Vector FGRDiffusionFlex::f1(int ib, int iFace, const Vector& h,
                        const Array3d<SymmetricTensor>& sigmaBlk)
{
   SymmetricTensor
      sigma = (sigmaBlk(ib) + sigmaBlk(ib+faceNbrOffset_[iFace])) / 2.0;
   Vector S(0, 0, 0);
   switch (iFace)
   {
     case 0:
      S[0] = -h[1]*h[2];
      break;
     case 1:
      S[1] = -h[0]*h[2];
      break;
     case 2:
      S[0] = h[1]*h[2];
      break;
     case 3:
      S[1] = h[0]*h[2];
      break;
     case 4:
      S[2] = -h[0]*h[1];
      break;
     case 5:
      S[2] = h[0]*h[1];
      break;
      
     default:
      assert(false);
   }
   return sigma * S;
}


void FGRDiffusionFlex::printAllWeights(const Array3d<int>& tissue)
{
   for (unsigned ii=0; ii<localTuple_.size(); ++ii)
   {
      Tuple globalTuple = localGrid_.globalTuple(localTuple_[ii]);
      int xx = localTuple_[ii].x();
      int yy = localTuple_[ii].y();
      int zz = localTuple_[ii].z();

      printf("DiffusionWeight: %5d %5d %5d %4d"
             " %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e"
             " %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e"
             " %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e"
             " %20.12e\n",
             globalTuple.x(), globalTuple.y(), globalTuple.z(),
             tissue(xx, yy, zz),
             weight_(xx, yy, zz).A[0],
             weight_(xx, yy, zz).A[1],
             weight_(xx, yy, zz).A[2],
             weight_(xx, yy, zz).A[3],
             weight_(xx, yy, zz).A[4],
             weight_(xx, yy, zz).A[5],
             weight_(xx, yy, zz).A[6],
             weight_(xx, yy, zz).A[7],
             weight_(xx, yy, zz).A[8],
             weight_(xx, yy, zz).A[9],
             weight_(xx, yy, zz).A[10],
             weight_(xx, yy, zz).A[11],
             weight_(xx, yy, zz).A[12],
             weight_(xx, yy, zz).A[13],
             weight_(xx, yy, zz).A[14],
             weight_(xx, yy, zz).A[15],
             weight_(xx, yy, zz).A[16],
             weight_(xx, yy, zz).A[17],
             weight_(xx, yy, zz).A[18]);
   }
}

void FGRDiffusionFlex::reorder_Coeff()
{
  uint32_t idx=0,ll,ii,jj,kk;
  const uint32_t Nx2 = nx;
  const uint32_t Ny2 = ny;
  const uint32_t Nz2 = nz;

  diffCoefT2_.resize(Nx2,Ny2,Nz2*19,0.0);
  srand(21321);

  for(int ii=0;ii<nLocal_;ii++)
  {
    int xx = localTuple_[ii].x();
    int yy = localTuple_[ii].y();
    int zz = localTuple_[ii].z();
    int zz4 = (int)(zz/4);
    int zp = zz + 1;
    int zm = zz - 1;
    int zp4 = (int)(zp/4);
    int zm4 = (int)(zm/4);

    assert(zz > 0);
    assert(zz < Nz2-1);

    diffCoefT2_(xx,yy,4*19*zp4 + 4* 0 + zp%4 ) =weight_(xx,yy,zz).A[ZZP];
    diffCoefT2_(xx,yy,4*19*zp4 + 4* 1 + zp%4 ) =weight_(xx,yy,zz).A[ZMP];
    diffCoefT2_(xx,yy,4*19*zp4 + 4* 2 + zp%4 ) =weight_(xx,yy,zz).A[MZP];
    diffCoefT2_(xx,yy,4*19*zp4 + 4* 3 + zp%4 ) =weight_(xx,yy,zz).A[ZPP];
    diffCoefT2_(xx,yy,4*19*zp4 + 4* 4 + zp%4 ) =weight_(xx,yy,zz).A[PZP];
                                 
    diffCoefT2_(xx,yy,4*19*zz4 + 4* 5 + zz%4 ) =weight_(xx,yy,zz).A[ZZZ];
    diffCoefT2_(xx,yy,4*19*zz4 + 4* 6 + zz%4 ) =weight_(xx,yy,zz).A[PMZ];
    diffCoefT2_(xx,yy,4*19*zz4 + 4* 7 + zz%4 ) =weight_(xx,yy,zz).A[MMZ];
    diffCoefT2_(xx,yy,4*19*zz4 + 4* 8 + zz%4 ) =weight_(xx,yy,zz).A[MPZ];
    diffCoefT2_(xx,yy,4*19*zz4 + 4* 9 + zz%4 ) =weight_(xx,yy,zz).A[PPZ];
    diffCoefT2_(xx,yy,4*19*zz4 + 4*10 + zz%4 ) =weight_(xx,yy,zz).A[ZMZ];
    diffCoefT2_(xx,yy,4*19*zz4 + 4*11 + zz%4 ) =weight_(xx,yy,zz).A[MZZ];
    diffCoefT2_(xx,yy,4*19*zz4 + 4*12 + zz%4 ) =weight_(xx,yy,zz).A[ZPZ];
    diffCoefT2_(xx,yy,4*19*zz4 + 4*13 + zz%4 ) =weight_(xx,yy,zz).A[PZZ];

    diffCoefT2_(xx,yy,4*19*zm4 + 4*14 + zm%4 ) =weight_(xx,yy,zz).A[ZZM];
    diffCoefT2_(xx,yy,4*19*zm4 + 4*15 + zm%4 ) =weight_(xx,yy,zz).A[ZMM];
    diffCoefT2_(xx,yy,4*19*zm4 + 4*16 + zm%4 ) =weight_(xx,yy,zz).A[MZM];
    diffCoefT2_(xx,yy,4*19*zm4 + 4*17 + zm%4 ) =weight_(xx,yy,zz).A[ZPM];
    diffCoefT2_(xx,yy,4*19*zm4 + 4*18 + zm%4 ) =weight_(xx,yy,zz).A[PZM];

  }
}

// simdiazed version
// 'start' cannot be in the middle. 'start' must point to (n,m,0). 
void
FGRDiffusionFlex::FGRDiff_simd_thread(const uint32_t b_quad,const int32_t e_quad, double* out0)
{
  int ii;
  double* VmM = VmBlock_.cBlock();

  vector4double B0,Sum1,B2,C0,C1,C2,Sum0,Sum2,Sum;
  vector4double my_x_y_z    ;
  vector4double my_xp1_y_z  ;
  vector4double my_x_yp1_z  ;
  vector4double my_xm1_y_z  ;
  vector4double my_x_ym1_z  ;
  vector4double my_xp1_yp1_z;
  vector4double my_xm1_yp1_z;
  vector4double my_xm1_ym1_z;
  vector4double my_xp1_ym1_z;
  vector4double my_zero_vec = vec_splats(0.0);

  double* out = out0 + b_quad*4;
  WeightType *simd_diff_ = &(diffCoefT3_[0]) + b_quad*19*4;
  uint16_t* vOffset = &(dOffset_[0]) + b_quad*9;

//    assert((uint64_t)phi_xm1_ym1_z%(4*sizeof(double)) == 0);
//    assert((uint64_t)phi_xm1_y_z  %(4*sizeof(double)) == 0);
//    assert((uint64_t)phi_xm1_yp1_z%(4*sizeof(double)) == 0);
//    assert((uint64_t)phi_x_ym1_z  %(4*sizeof(double)) == 0);
//    assert((uint64_t)phi_x_y_z    %(4*sizeof(double)) == 0);
//    assert((uint64_t)phi_x_yp1_z  %(4*sizeof(double)) == 0);
//    assert((uint64_t)phi_xp1_ym1_z%(4*sizeof(double)) == 0);
//    assert((uint64_t)phi_xp1_y_z  %(4*sizeof(double)) == 0);
//    assert((uint64_t)phi_xp1_yp1_z%(4*sizeof(double)) == 0);
//    assert((uint64_t)out          %(4*sizeof(double)) == 0);
 
  #define load_my_vectors  \
      my_x_y_z      =vec_ld(*(vOffset++), VmM);\
      my_x_ym1_z    =vec_ld(*(vOffset++), VmM);\
      my_x_yp1_z    =vec_ld(*(vOffset++), VmM);\
      my_xm1_y_z    =vec_ld(*(vOffset++), VmM);\
      my_xm1_ym1_z  =vec_ld(*(vOffset++), VmM);\
      my_xm1_yp1_z  =vec_ld(*(vOffset++), VmM);\
      my_xp1_y_z    =vec_ld(*(vOffset++), VmM);\
      my_xp1_ym1_z  =vec_ld(*(vOffset++), VmM);\
      my_xp1_yp1_z  =vec_ld(*(vOffset++), VmM);
 
   #define calc_zp(x) \
       x =  vec_madd(vec_ld(4*4*WTSZ,simd_diff_)  , my_xp1_y_z, \
            vec_madd(vec_ld(4*3*WTSZ,simd_diff_)  , my_x_yp1_z, \
            vec_madd(vec_ld(4*2*WTSZ,simd_diff_)  , my_xm1_y_z, \
            vec_madd(vec_ld(4*1*WTSZ,simd_diff_)  , my_x_ym1_z, \
            vec_mul( vec_ld(4*0*WTSZ,simd_diff_)  , my_x_y_z)))));
 
   #define calc_zz(x) \
       x = vec_madd( vec_ld(4*13*WTSZ,simd_diff_)  , my_xp1_y_z, \
           vec_madd( vec_ld(4*12*WTSZ,simd_diff_)  , my_x_yp1_z, \
           vec_madd( vec_ld(4*11*WTSZ,simd_diff_)  , my_xm1_y_z, \
           vec_madd( vec_ld(4*10*WTSZ,simd_diff_)  , my_x_ym1_z, \
           vec_madd( vec_ld(4*9*WTSZ,simd_diff_)  , my_xp1_yp1_z, \
           vec_madd( vec_ld(4*8*WTSZ,simd_diff_)  , my_xm1_yp1_z, \
           vec_madd( vec_ld(4*7*WTSZ,simd_diff_)  , my_xm1_ym1_z, \
           vec_madd( vec_ld(4*6*WTSZ,simd_diff_)  , my_xp1_ym1_z, \
           vec_mul ( vec_ld(4*5*WTSZ,simd_diff_)  , my_x_y_z)))))))));
 
   #define calc_zm(x) \
       x = vec_madd( vec_ld(4*18*WTSZ,simd_diff_)  , my_xp1_y_z, \
           vec_madd( vec_ld(4*17*WTSZ,simd_diff_)  , my_x_yp1_z, \
           vec_madd( vec_ld(4*16*WTSZ,simd_diff_)  , my_xm1_y_z, \
           vec_madd( vec_ld(4*15*WTSZ,simd_diff_)  , my_x_ym1_z, \
           vec_mul ( vec_ld(4*14*WTSZ,simd_diff_)  , my_x_y_z)))));
 
    load_my_vectors;
 
    calc_zp(B0);
    calc_zz(Sum1);
    calc_zm(B2);
 
    Sum2 = vec_sldw(my_zero_vec,B2,3);
 
  for(int qIdx=b_quad;qIdx<e_quad;qIdx++)
  {
    simd_diff_ += 19*4;
    load_my_vectors;
 
    calc_zp(C0);
    Sum0 = vec_sldw(B0,C0,1);
    B0=C0;
 
    Sum = vec_add(vec_add(Sum2,Sum1),Sum0);
    vec_st(Sum,0,out);  //commit
    out+=4;
 
    calc_zz(Sum1);
    calc_zm(C2);
    Sum2 = vec_sldw(B2,C2,3);
    B2 = C2;
  }
}

// simdiazed version
void
FGRDiffusionFlex::stencil(const uint32_t b_quad,const int32_t e_quad, double* out0)
{
  int ii;
  double* VmM = VmBlock_.cBlock();

  vector4double B0,Sum1,B2,C0,C1,C2,Sum0,Sum2,Sum;
  vector4double my_x_y_z    ;
  vector4double my_xp1_y_z  ;
  vector4double my_x_yp1_z  ;
  vector4double my_xm1_y_z  ;
  vector4double my_x_ym1_z  ;
  vector4double my_xp1_yp1_z;
  vector4double my_xm1_yp1_z;
  vector4double my_xm1_ym1_z;
  vector4double my_xp1_ym1_z;
  vector4double my_zero_vec = vec_splats(0.0);

  double* out = out0 + b_quad*4;
  WeightType *simd_diff_ = &(diffCoefT3_[0]) + b_quad*19*4;
  uint16_t* vOffset = &(dOffset_[0]) + b_quad*9;

//    assert((uint64_t)phi_xm1_ym1_z%(4*sizeof(double)) == 0);
//    assert((uint64_t)phi_xm1_y_z  %(4*sizeof(double)) == 0);
//    assert((uint64_t)phi_xm1_yp1_z%(4*sizeof(double)) == 0);
//    assert((uint64_t)phi_x_ym1_z  %(4*sizeof(double)) == 0);
//    assert((uint64_t)phi_x_y_z    %(4*sizeof(double)) == 0);
//    assert((uint64_t)phi_x_yp1_z  %(4*sizeof(double)) == 0);
//    assert((uint64_t)phi_xp1_ym1_z%(4*sizeof(double)) == 0);
//    assert((uint64_t)phi_xp1_y_z  %(4*sizeof(double)) == 0);
//    assert((uint64_t)phi_xp1_yp1_z%(4*sizeof(double)) == 0);
//    assert((uint64_t)out          %(4*sizeof(double)) == 0);
 
  #define load_my_vectors  \
      my_x_y_z      =vec_ld(*(vOffset++), VmM);\
      my_x_ym1_z    =vec_ld(*(vOffset++), VmM);\
      my_x_yp1_z    =vec_ld(*(vOffset++), VmM);\
      my_xm1_y_z    =vec_ld(*(vOffset++), VmM);\
      my_xm1_ym1_z  =vec_ld(*(vOffset++), VmM);\
      my_xm1_yp1_z  =vec_ld(*(vOffset++), VmM);\
      my_xp1_y_z    =vec_ld(*(vOffset++), VmM);\
      my_xp1_ym1_z  =vec_ld(*(vOffset++), VmM);\
      my_xp1_yp1_z  =vec_ld(*(vOffset++), VmM);
 
   #define calc_zp(x) \
       x =  vec_madd(vec_ld(4*4*WTSZ,simd_diff_)  , my_xp1_y_z, \
            vec_madd(vec_ld(4*3*WTSZ,simd_diff_)  , my_x_yp1_z, \
            vec_madd(vec_ld(4*2*WTSZ,simd_diff_)  , my_xm1_y_z, \
            vec_madd(vec_ld(4*1*WTSZ,simd_diff_)  , my_x_ym1_z, \
            vec_mul( vec_ld(4*0*WTSZ,simd_diff_)  , my_x_y_z)))));
 
   #define calc_zz(x) \
       x = vec_madd( vec_ld(4*13*WTSZ,simd_diff_)  , my_xp1_y_z, \
           vec_madd( vec_ld(4*12*WTSZ,simd_diff_)  , my_x_yp1_z, \
           vec_madd( vec_ld(4*11*WTSZ,simd_diff_)  , my_xm1_y_z, \
           vec_madd( vec_ld(4*10*WTSZ,simd_diff_)  , my_x_ym1_z, \
           vec_madd( vec_ld(4*9*WTSZ,simd_diff_)  , my_xp1_yp1_z, \
           vec_madd( vec_ld(4*8*WTSZ,simd_diff_)  , my_xm1_yp1_z, \
           vec_madd( vec_ld(4*7*WTSZ,simd_diff_)  , my_xm1_ym1_z, \
           vec_madd( vec_ld(4*6*WTSZ,simd_diff_)  , my_xp1_ym1_z, \
           vec_mul ( vec_ld(4*5*WTSZ,simd_diff_)  , my_x_y_z)))))))));
 
   #define calc_zm(x) \
       x = vec_madd( vec_ld(4*18*WTSZ,simd_diff_)  , my_xp1_y_z, \
           vec_madd( vec_ld(4*17*WTSZ,simd_diff_)  , my_x_yp1_z, \
           vec_madd( vec_ld(4*16*WTSZ,simd_diff_)  , my_xm1_y_z, \
           vec_madd( vec_ld(4*15*WTSZ,simd_diff_)  , my_x_ym1_z, \
           vec_mul ( vec_ld(4*14*WTSZ,simd_diff_)  , my_x_y_z)))));
 
    load_my_vectors;
    calc_zm(B2);

    simd_diff_ += 19*4;
    load_my_vectors;
    calc_zp(B0);
    calc_zz(Sum1);
    calc_zm(C2);
 
    Sum2 = vec_sldw(B2,C2,3);
    B2 = C2;
 
  for(int qIdx=b_quad;qIdx<e_quad;qIdx++)
  {
    simd_diff_ += 19*4;
    load_my_vectors;
 
    calc_zp(C0);
    Sum0 = vec_sldw(B0,C0,1);
    B0=C0;
 
    Sum = vec_add(vec_add(Sum2,Sum1),Sum0);
    vec_st(Sum,0,out);  //commit
    out+=4;
 
    calc_zz(Sum1);
    calc_zm(C2);
    Sum2 = vec_sldw(B2,C2,3);
    B2 = C2;
  }
}


void FGRDiffusionFlex::buildOffset(const Anatomy& anatomy)
{
   /////////////////////////////////////////////////////////////////////////////
   //simdization over arbitrary shape
   //////////////////////////////////////////////////////////////////////////////
   Array3d<int> boxToIdx;  //local + remote
   Array3d<int> cellToIdx; //local only

   assert(nz%4 == 0);
   int qnz=nz/4;

   QuadBlock_.resize(nx,ny,qnz,0);
   boxToIdx.resize(nx,ny,qnz,-1);
   cellToIdx.resize(nx,ny,qnz,-1);

   //mark cellness in QuadBlock_
   for(int ii=0;ii<anatomy.nLocal();ii++)
   {
     Tuple localTuple = localTuple_[ii];
//     if(anatomy.cellType(ii) < 3) //if cell  
     {
       int x = localTuple.x();
       int y = localTuple.y();
       int z = localTuple.z();
       int qz = z/4;

       QuadBlock_(x,y,qz) = 1;
       if (qz>0) { if (z%4 == 0) QuadBlock_(x,y,qz-1) = 1; }
       if (qz<qnz-1) {if ((z+1)%4 == 0) QuadBlock_(x,y,qz+1) = 1; }
     }
   }

   //halo mark
   for(int ii=anatomy.nLocal();ii<anatomy.size();ii++)
   {
     Tuple localTuple = localTuple_[ii];
//     if(anatomy.cellType(ii) < 3) //if cell  
     {
       int x = localTuple.x();
       int y = localTuple.y();
       int z = localTuple.z();
       int qz = z/4;

       if (QuadBlock_(x,y,qz) == 0) QuadBlock_(x,y,qz)=2 ;
     }
   }

   //idx build
   int bIdx=0;
   int cIdx=0;
   for(int ii=0;ii<nx;ii++)
   for(int jj=0;jj<ny;jj++)
   for(int kk=0;kk<qnz;kk++)
   {
     if (QuadBlock_(ii,jj,kk) != 0) boxToIdx(ii,jj,kk) = bIdx++;
     else boxToIdx(ii,jj,kk) = -1;

     if (QuadBlock_(ii,jj,kk) == 1) cellToIdx(ii,jj,kk) = cIdx++;
   }
   nQuad_=bIdx;
   nCalc_=cIdx;

//   cout << "nQuad:" << nQuad_ << endl;
//   cout << "nCalc:" << nCalc_ << endl;


   dOffset_.resize((nCalc_+1)*9);
   diffCoefT3_.resize((nCalc_+2)*4*19);

   //generate stream
   //go through quad marked as cell
   int dIdx=0;
   int eIdx=0;

   dOffset_[dIdx++]=U16_MAX;
   dOffset_[dIdx++]=U16_MAX;
   dOffset_[dIdx++]=U16_MAX;
   dOffset_[dIdx++]=U16_MAX;
   dOffset_[dIdx++]=U16_MAX;
   dOffset_[dIdx++]=U16_MAX;
   dOffset_[dIdx++]=U16_MAX;
   dOffset_[dIdx++]=U16_MAX;
   dOffset_[dIdx++]=U16_MAX;

   for(int ll=0;ll<4*19;ll++)
   {
     diffCoefT3_[eIdx++] = 0;
   }

   for(int ii=0;ii<nx-2;ii++)
   for(int jj=0;jj<ny-2;jj++)
   for(int kk=0;kk<qnz;kk++)
   {
     if (QuadBlock_(ii+1,jj+1,kk) == 1)
     {
//        dOffset_[dIdx++]=boxToIdx(ii+1,jj+1,kk)*4*8;
//        dOffset_[dIdx++]=boxToIdx(ii+1,jj+0,kk)*4*8;
//        dOffset_[dIdx++]=boxToIdx(ii+1,jj+2,kk)*4*8;
//        dOffset_[dIdx++]=boxToIdx(ii+0,jj+1,kk)*4*8;
//        dOffset_[dIdx++]=boxToIdx(ii+0,jj+0,kk)*4*8;
//        dOffset_[dIdx++]=boxToIdx(ii+0,jj+2,kk)*4*8;
//        dOffset_[dIdx++]=boxToIdx(ii+2,jj+1,kk)*4*8;
//        dOffset_[dIdx++]=boxToIdx(ii+2,jj+0,kk)*4*8;
//        dOffset_[dIdx++]=boxToIdx(ii+2,jj+2,kk)*4*8;

        dOffset_[dIdx++]=boxToIdx(ii+1,jj+1,kk)<0 ? U16_MAX:boxToIdx(ii+1,jj+1,kk)*4*8;
        dOffset_[dIdx++]=boxToIdx(ii+1,jj+0,kk)<0 ? U16_MAX:boxToIdx(ii+1,jj+0,kk)*4*8;
        dOffset_[dIdx++]=boxToIdx(ii+1,jj+2,kk)<0 ? U16_MAX:boxToIdx(ii+1,jj+2,kk)*4*8;
        dOffset_[dIdx++]=boxToIdx(ii+0,jj+1,kk)<0 ? U16_MAX:boxToIdx(ii+0,jj+1,kk)*4*8;
        dOffset_[dIdx++]=boxToIdx(ii+0,jj+0,kk)<0 ? U16_MAX:boxToIdx(ii+0,jj+0,kk)*4*8;
        dOffset_[dIdx++]=boxToIdx(ii+0,jj+2,kk)<0 ? U16_MAX:boxToIdx(ii+0,jj+2,kk)*4*8;
        dOffset_[dIdx++]=boxToIdx(ii+2,jj+1,kk)<0 ? U16_MAX:boxToIdx(ii+2,jj+1,kk)*4*8;
        dOffset_[dIdx++]=boxToIdx(ii+2,jj+0,kk)<0 ? U16_MAX:boxToIdx(ii+2,jj+0,kk)*4*8;
        dOffset_[dIdx++]=boxToIdx(ii+2,jj+2,kk)<0 ? U16_MAX:boxToIdx(ii+2,jj+2,kk)*4*8;

        for(int ll=0;ll<4*19;ll++)
        {
          diffCoefT3_[eIdx++] = diffCoefT2_(ii+1,jj+1,kk*4*19 + ll);
        }
     }
   }

   for(int ll=0;ll<4*19;ll++) { diffCoefT3_[eIdx++] = 0.0; }

   assert(dIdx%9 == 0);
   assert(dIdx/9 == (nCalc_+1));

   for(int ll=1;ll<9;ll++)
   {
     int mode=0;
     int last;
     for(int ii=0;ii<=nCalc_;ii++)
     {
       int jj=nCalc_-ii;
       //find first not -1
       if(mode==0)
       {
         if ((dOffset_[jj*9+ll]!=U16_MAX))
         {
           last = dOffset_[jj*9+ll]; 
           for(int kk=jj;kk<nCalc_;kk++) dOffset_[kk*9+ll]=last;
           mode=1;
         }
       }
       else
       {
         if (dOffset_[jj*9+ll] == U16_MAX ) dOffset_[jj*9+ll] = last;
         else last=dOffset_[jj*9+ll];
       }
    }
  }

  //Integrate copy index gen
  //Halo copy index gen
  inIndex_.resize(anatomy.size());
  for(int ii=0;ii<anatomy.size();ii++)
  {
    Tuple localTuple = localTuple_[ii];
//    if(anatomy.cellType(ii) < 3) //if cell  
    {
      int x = localTuple.x();
      int y = localTuple.y();
      int z = localTuple.z();
      int qz = z/4;

      inIndex_[ii] = boxToIdx(x,y,qz)*4 + z%4; 
    }
  }

  outIndex_.resize(anatomy.nLocal());
  for(int ii=0;ii<anatomy.nLocal();ii++)
  {
    Tuple localTuple = localTuple_[ii];
//    if(anatomy.cellType(ii) < 3) //if cell  
    {
      int x = localTuple.x();
      int y = localTuple.y();
      int z = localTuple.z();
      int qz = z/4;

      assert( (cellToIdx(x,y,qz)+1) > 0);
      outIndex_[ii] = cellToIdx(x,y,qz)*4 + z%4; 
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
// test functions
////////////////////////////////////////////////////////////////////////////////
void FGRDiffusionFlex::calc_history(int slice)
{
    std::cout << std::endl;
    for(int ii=0;ii<nCalc_;ii++)
    {
      std::cout << VmBlock_(dOffset_[9*ii+0]/8);
      std::cout << VmBlock_(dOffset_[9*ii+0]/8+1);
      std::cout << VmBlock_(dOffset_[9*ii+0]/8+2);
      std::cout << VmBlock_(dOffset_[9*ii+0]/8+3);
      std::cout << " " ;
      std::cout << VmBlock_(dOffset_[9*ii+1]/8);
      std::cout << VmBlock_(dOffset_[9*ii+1]/8+1);
      std::cout << VmBlock_(dOffset_[9*ii+1]/8+2);
      std::cout << VmBlock_(dOffset_[9*ii+1]/8+3);
      std::cout << " " ;
      std::cout << diffCoefT3_[ii*4*19  + 4*5 + 0];
      std::cout << diffCoefT3_[ii*4*19  + 4*5 + 1];
      std::cout << diffCoefT3_[ii*4*19  + 4*5 + 2];
      std::cout << diffCoefT3_[ii*4*19  + 4*5 + 3];
      std::cout << " " ;
//      std::cout << out[ii*4+0];
//      std::cout << out[ii*4+1];
//      std::cout << out[ii*4+2];
//      std::cout << out[ii*4+3];
      std::cout << std::endl;
    }

}



void FGRDiffusionFlex::VmBlockSet()
{
  srand(1234);

  cout << "VmBlock size=" << VmBlock_.size() << endl;

  for(int ii=0;ii<VmBlock_.size();ii++)
  {
    VmBlock_(ii) = rand()%9 + 1;
  }
}


void FGRDiffusionFlex::dump_anatomy(const Anatomy& anatomy, int slice)
{

  for(int ii=0;ii<nx;ii++)
  for(int jj=0;jj<ny;jj++)
  for(int kk=0;kk<nz;kk++)
    VmBlock_(ii,jj,kk)=0;

  for(int ii=0;ii<anatomy.size();ii++)
  {
    Tuple localTuple = localTuple_[ii];
//    if(anatomy.cellType(ii) < 3) //if cell  
    {
      int x = localTuple.x();
      int y = localTuple.y();
      int z = localTuple.z();

      if(ii < nLocal_)
        VmBlock_(x,y,z) = 1;
      else
        VmBlock_(x,y,z) = (VmBlock_(x,y,z) == 0) ? 2:VmBlock_(x,y,z);
    }
  }

  srand(1234);
  //random assign weight
  for(int ii=0;ii<anatomy.nLocal();ii++)
  {
//    if(anatomy.cellType(ii) < 3) //if cell  
    {
      for(int jj=0;jj<19;jj++)
      {
        weight_(blockIndex_[ii]).A[jj] = 0;
//        if( jj == 5 )
          if ( VmBlock_(blockIndex_[ii] + offset_[jj] ) != 0 ) 
            weight_(blockIndex_[ii]).A[jj] = rand()%3+1;
      }
    }
  }
  
  dump_array3d(VmBlock_,slice);

}

void FGRDiffusionFlex::dump_quad_anatomy(const Anatomy& anatomy, int slice)
{

  for(int ii=0;ii<nx;ii++)
  for(int jj=0;jj<ny;jj++)
  for(int kk=0;kk<nz;kk++)
  {
    VmBlock_(ii,jj,kk)=QuadBlock_(ii,jj,kk/4);
  }

  dump_array3d(VmBlock_,slice);
}

void FGRDiffusionFlex::dump_array3d(Array3d<double>& Box,int slice)
{
  cout << endl;
  for(int ii=0;ii<nx;ii++)
  {
    for(int jj=0;jj<nz;jj++)
    {
      cout << Box(ii,slice,jj) << " " ;
    }
    cout << endl;
  }
}


void FGRDiffusionFlex::boundary_test(const Anatomy& anatomy)
{
  Array3d<int>                    anatomy_3d;
  srand(3245);

  anatomy_3d.resize(nx,ny,nz,0);
  //anatomy creation
  for(int ii=0;ii<anatomy.size();ii++)
  {
    anatomy_3d(blockIndex_[ii])=1;
  }

  //check if weight is zero for non-cell
  cout << "WeightSumTolerance=" << weightSumTolerance << endl;
  for(int ii=0;ii<nLocal_;ii++)
  {
    for(int jj=0;jj<19;jj++)
    {
      if(anatomy_3d(blockIndex_[ii]+offset_[jj]) == 0)
      {
        if(weight_(blockIndex_[ii]).A[jj] > weightSumTolerance)
        {
           cout << "Warning:imposing zero weight at " << ii << "," << jj <<  ". it was " << weight_(blockIndex_[ii]).A[jj] << endl;
           weight_(blockIndex_[ii]).A[jj]=0.0;
        }
      }
      else
      {
//           weight_(blockIndex_[ii]).A[jj]=0;
//           if(jj==1)  weight_(blockIndex_[ii]).A[1]=0.001;
           ;
           
      }
    }
  }
}

void FGRDiffusionFlex::basic_diffusion(Array3d<FGRUtils::DiffWeight>& coef, Array3d<double> &inV, Array3d<double> &outV)
{
  for (unsigned ii=0; ii<nLocal_; ++ii)
  {
    int idx = blockIndex_[ii];
    double sum=0;
    for(int kk=0;kk<19;kk++)
    {
      sum += inV(idx + offset_[kk]) * coef(idx).A[kk];
//      if ( ii==2) cout << kk << ":" << inV(idx + offset_[kk]) << " " << weight_(idx).A[kk] <<":";
    }
//     if(ii==2) cout << sum << endl;
    outV(idx) = sum;
  }
  cout << endl;
}
