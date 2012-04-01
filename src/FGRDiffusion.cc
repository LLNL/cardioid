#include "FGRDiffusion.hh"
#include <cassert>
#include "DiffusionUtils.hh"
#include "Anatomy.hh"
#include "Vector.hh"
#include "fastBarrier.hh"
#include <algorithm>
#include <cstdio>
#include "ThreadServer.hh"
#ifdef TIMING
#include "PerformanceTimers.hh"
using namespace PerformanceTimers;
#endif

//#define check_same
using namespace std;
using namespace FGRUtils;


FGRDiffusion::FGRDiffusion(const FGRDiffusionParms& parms,
                           const Anatomy& anatomy,
                           const ThreadTeam& threadInfo)
: localGrid_(DiffusionUtils::findBoundingBox_simd(anatomy)),
  threadInfo_(threadInfo),
  diffusionScale_(parms.diffusionScale_)
{

   unsigned nx = localGrid_.nx();
   unsigned ny = localGrid_.ny();
   unsigned nz = localGrid_.nz();

   // This is a test
   for (unsigned ii=0; ii<anatomy.size(); ++ii)
   {
      Tuple globalTuple = anatomy.globalTuple(ii);
      Tuple ll = localGrid_.localTuple(globalTuple);
      assert(ll.x() >= 0 && ll.y() >= 0 && ll.z() >= 0);
      assert(ll.x() < nx && ll.y() < ny && ll.z() < nz);
   }
   // This has been a test

   int chunkSize = anatomy.nLocal() / threadInfo.nThreads();
   int leftOver  = anatomy.nLocal() % threadInfo.nThreads();
   threadOffset_.resize(threadInfo.nThreads()+1);
   threadOffset_[0] = 0;
   for (int ii=0; ii<threadInfo.nThreads(); ++ii)
   {
      threadOffset_[ii+1] = threadOffset_[ii] + chunkSize;
      if (ii < leftOver)
         ++threadOffset_[ii+1];
   }
   assert(threadOffset_[threadInfo_.nThreads()] == anatomy.nLocal());

   barrierHandle_.resize(threadInfo.nThreads());
   fgrBarrier_ = L2_BarrierWithSync_InitShared();
   #pragma omp parallel
   {
      int tid = threadInfo.teamRank();
      if (tid >= 0)
         L2_BarrierWithSync_InitInThread(fgrBarrier_, &barrierHandle_[tid]);
   }         
   
   weight_.resize(nx, ny, nz);
   VmBlock_.resize(nx, ny, nz);

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
   reorder_Coeff();

   tmp_dVm.resize(nx,ny,nz + (nz%4==0 ? 0:4-(nz%4)));

   //simd thread offsets
   threadOffsetSimd_.resize(threadInfo.nThreads()+1);
   int tid=0;
   for (int ii=0; ii<threadInfo.nThreads(); ++ii) threadOffsetSimd_[ii+1]=0;
   for (int ii=1; ii<(nx-1); ++ii)  //1 ... nx-2
   {
      threadOffsetSimd_[tid+1]++;
      tid++; tid=tid % threadInfo.nThreads();
   }
   //printf("thread:1 ");
   threadOffsetSimd_[0]=1;
   for (int ii=0; ii<threadInfo.nThreads(); ++ii)
   {
      //printf("%d ",threadOffsetSimd_[ii+1]);
      threadOffsetSimd_[ii+1] += threadOffsetSimd_[ii];
   }
   //printf("\n");
   assert(nx-1 == threadOffsetSimd_[threadInfo.nThreads()] );
}




/** threaded simd version */
void FGRDiffusion::calc(const vector<double>& Vm, vector<double>& dVm, double *recv_buf_, int nLocal)
{
   int tid = threadInfo_.teamRank();
#ifdef TIMING
   profileStart(FGR_Array2MatrixTimer);
#endif
   if ( tid==0 ) updateVoltageBlock(Vm, recv_buf_, nLocal);
#ifdef TIMING
   profileStop(FGR_Array2MatrixTimer);
   profileStart(FGR_BarrierTimer);
#endif
   L2_BarrierWithSync_Barrier(fgrBarrier_, &barrierHandle_[tid], threadInfo_.nThreads());
#ifdef TIMING
   profileStop(FGR_BarrierTimer);
   profileStart(FGR_AlignCopyTimer);
#endif

   Array3d<double> *VmTmp = &(VmBlock_);

   //make sure z is multiple of 4
   if(VmBlock_.nz()%4 != 0)
   {
     int nz_4 = VmBlock_.nz() + 4-VmBlock_.nz()%4;
     VmTmp = new Array3d<double>;
     VmTmp->resize(VmBlock_.nx(),VmBlock_.ny(),nz_4);
     for(int ii=0;ii<VmBlock_.nx();ii++)
     for(int jj=0;jj<VmBlock_.ny();jj++)
     for(int kk=0;kk<nz_4;kk++)
     {
       (*VmTmp)(ii,jj,kk) = VmBlock_(ii,jj,kk);
     }
   }

#ifdef TIMING
   profileStop(FGR_AlignCopyTimer);
#endif
   
#ifdef TIMING
   profileStart(FGR_StencilTimer);
#endif
   //uint32_t begin = VmTmp->tupleToIndex(threadOffsetSimd_[tid],1,0);
   //uint32_t end = VmTmp->tupleToIndex(threadOffsetSimd_[tid+1]-1,VmTmp->ny()-2,VmTmp->nz());
   //   printf("simd version:%d-%d\n",begin,end);
   //if (threadOffsetSimd_[tid] < threadOffsetSimd_[tid+1] )
   //   FGRDiff_simd_thread(begin,end-begin,VmTmp,tmp_dVm.cBlock());
   if (threadOffsetSimd_[tid] < threadOffsetSimd_[tid+1] )
      FGRDiff_simd_thread(threadOffsetSimd_[tid] ,threadOffsetSimd_[tid+1],VmTmp,tmp_dVm.cBlock());

#ifdef TIMING
   profileStop(FGR_StencilTimer);
   profileStart(FGR_Barrier2Timer);
#endif
   L2_BarrierWithSync_Barrier(fgrBarrier_, &barrierHandle_[tid], threadInfo_.nThreads());
#ifdef TIMING
   profileStop(FGR_Barrier2Timer);
   profileStart(FGR_Matrix2ArrayTimer);
#endif
   if(VmBlock_.nz()%4 != 0) delete VmTmp;

   int begin = threadOffset_[tid];
   int end   = threadOffset_[tid+1];
   for (int ii=begin; ii<end; ++ii)
   {
      dVm[ii] = tmp_dVm(localTuple_[ii].x(),localTuple_[ii].y(),localTuple_[ii].z());
      dVm[ii] *= diffusionScale_;
   }
#ifdef TIMING
   profileStop(FGR_Matrix2ArrayTimer);
#endif
}

/** We're building the localTuple array only for local cells.  We can't
 * do stencil operations on remote particles so we shouldn't need
 * tuples.  We can use block indices instead.
 */
void FGRDiffusion::buildTupleArray(const Anatomy& anatomy)
{
   localTuple_.resize(anatomy.nLocal(), Tuple(0,0,0));
   for (unsigned ii=0; ii<anatomy.nLocal(); ++ii)
   {
      Tuple globalTuple = anatomy.globalTuple(ii);
      localTuple_[ii] = localGrid_.localTuple(globalTuple);
      assert(localTuple_[ii].x() > 0);
      assert(localTuple_[ii].y() > 0);
      assert(localTuple_[ii].z() > 0);
      assert(localTuple_[ii].x() < localGrid_.nx()-1);
      assert(localTuple_[ii].y() < localGrid_.ny()-1);
      assert(localTuple_[ii].z() < localGrid_.nz()-1);
   }
}

void FGRDiffusion::buildBlockIndex(const Anatomy& anatomy)
{
   blockIndex_.resize(anatomy.size());
   for (unsigned ii=0; ii<anatomy.size(); ++ii)
   {
      Tuple globalTuple = anatomy.globalTuple(ii);
      Tuple ll = localGrid_.localTuple(globalTuple);
      blockIndex_[ii] = VmBlock_.tupleToIndex(ll.x(), ll.y(), ll.z());
   }
}

void FGRDiffusion::precomputeCoefficients(const Anatomy& anatomy)
{
   unsigned nx = localGrid_.nx();
   unsigned ny = localGrid_.ny();
   unsigned nz = localGrid_.nz();
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
      assert(abs(sum) < 1e-14);
      
   }
//   printAllWeights(tissueBlk);
}


void FGRDiffusion::updateVoltageBlock(const vector<double>& Vm, double *recv_buf, int nLocal)
{
   for (unsigned ii=0; ii<nLocal; ++ii)
   {
      int index = blockIndex_[ii];
      VmBlock_(index) = Vm[ii];
   #ifdef check_same
      VmBlock_(index) = rand();
   #endif
   }
   assert(nLocal <= Vm.size());
   for (unsigned ii=nLocal; ii<Vm.size(); ++ii)
   {
      int index = blockIndex_[ii];
      VmBlock_(index) = recv_buf[ii-nLocal];
   #ifdef check_same
      VmBlock_(index) = rand();
   #endif
   }
}


void FGRDiffusion::mkTissueArray(
   const Array3d<int>& tissueBlk, int ib, int* tissue)
{
   for (unsigned ii=0; ii<19; ++ii)
      tissue[ii] = tissueBlk(ib + offset_[ii]);
}


Vector FGRDiffusion::f1(int ib, int iFace, const Vector& h,
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


void FGRDiffusion::printAllWeights(const Array3d<int>& tissue)
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

void FGRDiffusion::reorder_Coeff()
{
  uint32_t idx=0,ll,ii,jj,kk;
  const uint32_t Nx2 = localGrid_.nx();
  const uint32_t Ny2 = localGrid_.ny();
  const uint32_t Nz2 = localGrid_.nz() + ((localGrid_.nz()%4)==0 ? 0 : (4-localGrid_.nz()%4));

  diffCoefT2_.resize(Nx2,Ny2,Nz2*19,0.0);
  srand(21321);

  const int n_cell = localTuple_.size() ;
  for(int ii=0;ii<n_cell;ii++)
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
FGRDiffusion::FGRDiff_simd_thread(const uint32_t bx,const int32_t ex, Array3d<double>* VmTmp, double* out0)
{
  int ii;

  const unsigned Nx2 = VmTmp->nx();
  const unsigned Ny2 = VmTmp->ny();
  const unsigned Nz2 = VmTmp->nz();

  double* VmM = VmTmp->cBlock();
  double* out;
  double* diffC = diffCoefT2_.cBlock();

  int xm1ym1z_ = ((-1) *Ny2 + (-1)) * Nz2;
  int xm1yz_ =   ((-1) *Ny2 + ( 0)) * Nz2;
  int xm1yp1z_ = ((-1) *Ny2 + (+1)) * Nz2;
  int xym1z_ =   ((0 ) *Ny2 + (-1)) * Nz2;
  int xyp1z_ =   ((0 ) *Ny2 + (+1)) * Nz2;
  int xp1ym1z_ = ((+1) *Ny2 + (-1)) * Nz2;
  int xp1yz_ =   ((+1) *Ny2 + ( 0)) * Nz2;
  int xp1yp1z_ = ((+1) *Ny2 + (+1)) * Nz2;

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

  for(int xx=bx;xx<ex;xx++)
  {
    int start = VmTmp->tupleToIndex(xx,1,0);
    int end   = VmTmp->tupleToIndex(xx,Ny2-2,Nz2);
    int chunk_size = end - start;

    double* phi_xm1_ym1_z = VmM + start + xm1ym1z_;
    double* phi_xm1_y_z   = VmM + start + xm1yz_;
    double* phi_xm1_yp1_z = VmM + start + xm1yp1z_;
    double* phi_x_ym1_z   = VmM + start + xym1z_;
    double* phi_x_y_z     = VmM + start ;
    double* phi_x_yp1_z   = VmM + start + xyp1z_;
    double* phi_xp1_ym1_z = VmM + start + xp1ym1z_;
    double* phi_xp1_y_z   = VmM + start + xp1yz_;
    double* phi_xp1_yp1_z = VmM + start + xp1yp1z_;
 
    out = out0 + start;
 
    double *simd_diff_ = diffC + start*19; 
 
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
        my_x_y_z      =vec_ld(0,      phi_x_y_z   );\
        my_xp1_y_z    =vec_ld(0,      phi_xp1_y_z );\
        my_x_yp1_z    =vec_ld(0,      phi_x_yp1_z );\
        my_xm1_y_z    =vec_ld(0,      phi_xm1_y_z );\
        my_x_ym1_z    =vec_ld(0,      phi_x_ym1_z );\
        my_xp1_yp1_z  =vec_ld(0,    phi_xp1_yp1_z );\
        my_xm1_yp1_z  =vec_ld(0,    phi_xm1_yp1_z );\
        my_xm1_ym1_z  =vec_ld(0,    phi_xm1_ym1_z );\
        my_xp1_ym1_z  =vec_ld(0,    phi_xp1_ym1_z );
 
    #define calc_zp(x) \
        x =  vec_madd(vec_ld(4*4*8,simd_diff_)  , my_xp1_y_z, \
             vec_madd(vec_ld(4*3*8,simd_diff_)  , my_x_yp1_z, \
             vec_madd(vec_ld(4*2*8,simd_diff_)  , my_xm1_y_z, \
             vec_madd(vec_ld(4*1*8,simd_diff_)  , my_x_ym1_z, \
             vec_mul( vec_ld(4*0*8,simd_diff_)  , my_x_y_z)))));
 
    #define calc_zz(x) \
        x = vec_madd( vec_ld(4*13*8,simd_diff_)  , my_xp1_y_z, \
            vec_madd( vec_ld(4*12*8,simd_diff_)  , my_x_yp1_z, \
            vec_madd( vec_ld(4*11*8,simd_diff_)  , my_xm1_y_z, \
            vec_madd( vec_ld(4*10*8,simd_diff_)  , my_x_ym1_z, \
            vec_madd( vec_ld(4*9*8,simd_diff_)  , my_xp1_yp1_z, \
            vec_madd( vec_ld(4*8*8,simd_diff_)  , my_xm1_yp1_z, \
            vec_madd( vec_ld(4*7*8,simd_diff_)  , my_xm1_ym1_z, \
            vec_madd( vec_ld(4*6*8,simd_diff_)  , my_xp1_ym1_z, \
            vec_mul ( vec_ld(4*5*8,simd_diff_)  , my_x_y_z)))))))));
 
    #define calc_zm(x) \
        x = vec_madd( vec_ld(4*18*8,simd_diff_)  , my_xp1_y_z, \
            vec_madd( vec_ld(4*17*8,simd_diff_)  , my_x_yp1_z, \
            vec_madd( vec_ld(4*16*8,simd_diff_)  , my_xm1_y_z, \
            vec_madd( vec_ld(4*15*8,simd_diff_)  , my_x_ym1_z, \
            vec_mul ( vec_ld(4*14*8,simd_diff_)  , my_x_y_z)))));
 
    #define shift_pointers \
            phi_xp1_y_z   +=4; \
            phi_x_yp1_z   +=4; \
            phi_xm1_y_z   +=4; \
            phi_x_ym1_z   +=4; \
            phi_xp1_yp1_z +=4; \
            phi_xm1_yp1_z +=4; \
            phi_xm1_ym1_z +=4; \
            phi_xp1_ym1_z +=4; \
            phi_x_y_z     +=4;
 
    load_my_vectors;
 
    calc_zp(B0);
    calc_zz(Sum1);
    calc_zm(B2);
 
    Sum2 = vec_sldw(my_zero_vec,B2,3);
 
    for (ii=0;ii<chunk_size;ii+=4)
    {
      simd_diff_ += 19*4;
      shift_pointers; //note that these pointers one step further
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
}

// simdiazed version
// 'start' cannot be in the middle. 'start' must point to (n,0,0). 
void
FGRDiffusion::FGRDiff_simd(const uint32_t start,const int32_t chunk_size, Array3d<double>* VmTmp, double* out)
{
  int ii;

  const unsigned Nx2 = VmTmp->nx();
  const unsigned Ny2 = VmTmp->ny();
  const unsigned Nz2 = VmTmp->nz();

  double* VmM = VmTmp->cBlock();

  int xm1ym1z_ = ((-1) *Ny2 + (-1)) * Nz2;
  int xm1yz_ =   ((-1) *Ny2 + ( 0)) * Nz2;
  int xm1yp1z_ = ((-1) *Ny2 + (+1)) * Nz2;
  int xym1z_ =   ((0 ) *Ny2 + (-1)) * Nz2;
  int xyp1z_ =   ((0 ) *Ny2 + (+1)) * Nz2;
  int xp1ym1z_ = ((+1) *Ny2 + (-1)) * Nz2;
  int xp1yz_ =   ((+1) *Ny2 + ( 0)) * Nz2;
  int xp1yp1z_ = ((+1) *Ny2 + (+1)) * Nz2;

  double* phi_xm1_ym1_z = VmM + start + xm1ym1z_;
  double* phi_xm1_y_z   = VmM + start + xm1yz_;
  double* phi_xm1_yp1_z = VmM + start + xm1yp1z_;
  double* phi_x_ym1_z   = VmM + start + xym1z_;
  double* phi_x_y_z     = VmM + start ;
  double* phi_x_yp1_z   = VmM + start + xyp1z_;
  double* phi_xp1_ym1_z = VmM + start + xp1ym1z_;
  double* phi_xp1_y_z   = VmM + start + xp1yz_;
  double* phi_xp1_yp1_z = VmM + start + xp1yp1z_;

  out += start;

  assert((uint64_t)phi_xm1_ym1_z%(4*sizeof(double)) == 0);
  assert((uint64_t)phi_xm1_y_z  %(4*sizeof(double)) == 0);
  assert((uint64_t)phi_xm1_yp1_z%(4*sizeof(double)) == 0);
  assert((uint64_t)phi_x_ym1_z  %(4*sizeof(double)) == 0);
  assert((uint64_t)phi_x_y_z    %(4*sizeof(double)) == 0);
  assert((uint64_t)phi_x_yp1_z  %(4*sizeof(double)) == 0);
  assert((uint64_t)phi_xp1_ym1_z%(4*sizeof(double)) == 0);
  assert((uint64_t)phi_xp1_y_z  %(4*sizeof(double)) == 0);
  assert((uint64_t)phi_xp1_yp1_z%(4*sizeof(double)) == 0);
  assert((uint64_t)out          %(4*sizeof(double)) == 0);

  //initial
  vector4double B0,B1,B2,C0,C1,C2,Sum0,Sum2,Sum;

  vector4double my_x_y_z      =vec_ld(0,      phi_x_y_z   );
  vector4double my_xp1_y_z    =vec_ld(0,      phi_xp1_y_z );
  vector4double my_x_yp1_z    =vec_ld(0,      phi_x_yp1_z );
  vector4double my_xm1_y_z    =vec_ld(0,      phi_xm1_y_z );
  vector4double my_x_ym1_z    =vec_ld(0,      phi_x_ym1_z );
  vector4double my_xp1_yp1_z  =vec_ld(0,    phi_xp1_yp1_z );
  vector4double my_xm1_yp1_z  =vec_ld(0,    phi_xm1_yp1_z );
  vector4double my_xm1_ym1_z  =vec_ld(0,    phi_xm1_ym1_z );
  vector4double my_xp1_ym1_z  =vec_ld(0,    phi_xp1_ym1_z );
  vector4double my_zero_vec = vec_splats(0.0);
 
  double *simd_diff_ = diffCoefT2_.cBlock() + start*19;

  B0 =  vec_madd(vec_ld(4*4*8,simd_diff_)  , my_xp1_y_z,
        vec_madd(vec_ld(3*4*8,simd_diff_)  , my_x_yp1_z,
        vec_madd(vec_ld(2*4*8,simd_diff_)  , my_xm1_y_z,
        vec_madd(vec_ld(1*4*8,simd_diff_)  , my_x_ym1_z,
        vec_mul( vec_ld(0*4*8,simd_diff_)  , my_x_y_z)))));

  B1 = vec_madd( vec_ld(13*4*8,simd_diff_)  , my_xp1_y_z,
       vec_madd( vec_ld(12*4*8,simd_diff_)  , my_x_yp1_z,
       vec_madd( vec_ld(11*4*8,simd_diff_)  , my_xm1_y_z,
       vec_madd( vec_ld(10*4*8,simd_diff_)  , my_x_ym1_z,
       vec_madd( vec_ld( 9*4*8,simd_diff_)  , my_xp1_yp1_z,
       vec_madd( vec_ld( 8*4*8,simd_diff_)  , my_xm1_yp1_z,
       vec_madd( vec_ld( 7*4*8,simd_diff_)  , my_xm1_ym1_z,
       vec_madd( vec_ld( 6*4*8,simd_diff_)  , my_xp1_ym1_z,
       vec_mul ( vec_ld( 5*4*8,simd_diff_)  , my_x_y_z)))))))));

  B2 = vec_madd( vec_ld(18*4*8,simd_diff_)  , my_xp1_y_z,
       vec_madd( vec_ld(17*4*8,simd_diff_)  , my_x_yp1_z,
       vec_madd( vec_ld(16*4*8,simd_diff_)  , my_xm1_y_z,
       vec_madd( vec_ld(15*4*8,simd_diff_)  , my_x_ym1_z,
       vec_mul ( vec_ld(14*4*8,simd_diff_)  , my_x_y_z)))));

  Sum2 = vec_sldw(my_zero_vec,B2,3);

  for(ii=0;ii<chunk_size;ii+=4)  //start must be that os x chunk.
  {
    simd_diff_ += 19*4;
    phi_xp1_y_z   +=4;
    phi_x_yp1_z   +=4;
    phi_xm1_y_z   +=4;
    phi_x_ym1_z   +=4;
    phi_xp1_yp1_z +=4;
    phi_xm1_yp1_z +=4;
    phi_xm1_ym1_z +=4;
    phi_xp1_ym1_z +=4;
    phi_x_y_z     +=4;

    my_x_y_z      =vec_ld(0,      phi_x_y_z   );
    my_xp1_y_z    =vec_ld(0,      phi_xp1_y_z );
    my_x_yp1_z    =vec_ld(0,      phi_x_yp1_z );
    my_xm1_y_z    =vec_ld(0,      phi_xm1_y_z );
    my_x_ym1_z    =vec_ld(0,      phi_x_ym1_z );
    my_xp1_yp1_z  =vec_ld(0,    phi_xp1_yp1_z );
    my_xm1_yp1_z  =vec_ld(0,    phi_xm1_yp1_z );
    my_xm1_ym1_z  =vec_ld(0,    phi_xm1_ym1_z );
    my_xp1_ym1_z  =vec_ld(0,    phi_xp1_ym1_z );


  C0 =  vec_madd(vec_ld(4*4*8,simd_diff_)  , my_xp1_y_z,
        vec_madd(vec_ld(3*4*8,simd_diff_)  , my_x_yp1_z,
        vec_madd(vec_ld(2*4*8,simd_diff_)  , my_xm1_y_z,
        vec_madd(vec_ld(1*4*8,simd_diff_)  , my_x_ym1_z,
        vec_mul( vec_ld(0*4*8,simd_diff_)  , my_x_y_z)))));

    Sum0 = vec_sldw(B0,C0,1);
    B0 = C0;

    // commit the previous
    Sum = vec_add(vec_add(Sum0,B1),Sum2);

    vec_st(Sum,0,out);

  B1 = vec_madd( vec_ld(13*4*8,simd_diff_)  , my_xp1_y_z,
       vec_madd( vec_ld(12*4*8,simd_diff_)  , my_x_yp1_z,
       vec_madd( vec_ld(11*4*8,simd_diff_)  , my_xm1_y_z,
       vec_madd( vec_ld(10*4*8,simd_diff_)  , my_x_ym1_z,
       vec_madd( vec_ld( 9*4*8,simd_diff_)  , my_xp1_yp1_z,
       vec_madd( vec_ld( 8*4*8,simd_diff_)  , my_xm1_yp1_z,
       vec_madd( vec_ld( 7*4*8,simd_diff_)  , my_xm1_ym1_z,
       vec_madd( vec_ld( 6*4*8,simd_diff_)  , my_xp1_ym1_z,
       vec_mul ( vec_ld( 5*4*8,simd_diff_)  , my_x_y_z)))))))));

  C2 = vec_madd( vec_ld(18*4*8,simd_diff_)  , my_xp1_y_z,
       vec_madd( vec_ld(17*4*8,simd_diff_)  , my_x_yp1_z,
       vec_madd( vec_ld(16*4*8,simd_diff_)  , my_xm1_y_z,
       vec_madd( vec_ld(15*4*8,simd_diff_)  , my_x_ym1_z,
       vec_mul ( vec_ld(14*4*8,simd_diff_)  , my_x_y_z)))));

    Sum2 = vec_sldw(B2,C2,3);
    B2 = C2;

   out +=4;

  }
}

