#include "FGRDiffusionOverlap.hh"
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

//#define DIFF_TEST

FGRDiffusionOverlap::FGRDiffusionOverlap(const FGRDiffusionParms& parms,
                           const Anatomy& anatomy,
                           const ThreadTeam& threadInfo,
                           const ThreadTeam& reactionThreadInfo)
: nLocal_(anatomy.nLocal()),nRemote_(anatomy.nRemote()),
  localGrid_(DiffusionUtils::findBoundingBox_simd(anatomy, parms.printBBox_)),
  threadInfo_(threadInfo),
  reactionThreadInfo_(reactionThreadInfo),
  diffusionScale_(parms.diffusionScale_)
{

   unsigned nx = localGrid_.nx();
   unsigned ny = localGrid_.ny();
   unsigned nz = localGrid_.nz();




   cout << "the bounding box size is"<<nx<<"x"<<ny<<"x"<<nz<<endl;
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

  
   //simd thread offsets
   int work[3]={ny,nx,nx};
   for(int ii=0;ii<3;ii++)
   {
     int Twork=work[ii]-2;
     int chunkSize = Twork / threadInfo.nThreads();
     int leftOver = Twork % threadInfo.nThreads();
     threadOffset2D_[ii].resize(threadInfo.nThreads()+1);
     threadOffset2D_[ii][0]=1;
     //cout << "2D trhead offset:";
     for (int jj=0; jj<threadInfo.nThreads(); ++jj)
     {
        threadOffset2D_[ii][jj+1] = threadOffset2D_[ii][jj] + chunkSize;
        if (jj < leftOver)
           ++threadOffset2D_[ii][jj+1];
        //cout <<ii<<":"<<jj<<":"<<threadOffset2D_[ii][jj+1]<<":";
     }
     //cout << endl;
     assert(work[ii]-1 == threadOffset2D_[ii][threadInfo.nThreads()] );
   }
     
   int chunkSize = (nx-2) / threadInfo.nThreads();
   int leftOver = (nx-2) % threadInfo.nThreads();
   threadOffsetSimd_.resize(threadInfo.nThreads()+1);
   threadOffsetSimd_[0]=1;
   for (int ii=0; ii<threadInfo.nThreads(); ++ii)
   {
      threadOffsetSimd_[ii+1] = threadOffsetSimd_[ii] + chunkSize;
      if (ii < leftOver)
         ++threadOffsetSimd_[ii+1];
   }
   assert(nx-1 == threadOffsetSimd_[threadInfo.nThreads()] );

   barrierHandle_.resize(threadInfo.nThreads());
   fgrBarrier_ = L2_BarrierWithSync_InitShared();
   #pragma omp parallel
   {
      int tid = threadInfo.teamRank();
      if (tid >= 0)
         L2_BarrierWithSync_InitInThread(fgrBarrier_, &barrierHandle_[tid]);
   }         
   
   weight_.resize(nx, ny, nz);
   A0_.resize(nx, ny, nz, 0.0);
   assert(nz%4 == 0);

   VmBlock_.resize(nx,ny,nz,0.0);
   VmSlab_[0].resize(1,ny,nz,0.0);
   VmSlab_[1].resize(1,ny,nz,0.0);
   VmSlab_[2].resize(1,nx,nz,0.0);
   VmSlab_[3].resize(1,nx,nz,0.0);
   VmSlab_[4].resize(1,nx,((int)(ny+3)/4)*4,0.0);
   VmSlab_[5].resize(1,nx,((int)(ny+3)/4)*4,0.0);

   dVmBlock_.resize(nx,ny,nz,0.0);

   buildTupleArray(anatomy);
   buildBlockIndex(anatomy);

   dVmBlock2Doffset_[0] =  ny*nz;   //x=1 plane;
   dVmBlock2Doffset_[1] =  (nx-2)*ny*nz;   //x=nx-2 plane;
   dVmBlock2Doffset_[2] =  nz  ;   //(0,1,0)//y=1 plane;
   dVmBlock2Doffset_[3] =  (ny-2)*nz;   //y=ny-2 plane;
   dVmBlock2Doffset_[4] =  BBzb+1;   //z=1 plane;
   dVmBlock2Doffset_[5] =  BBze-1;   //z=nz-1 plane;

   dVmBlock2Djump_[0] = dVmBlock2Djump_[1] = nz;
   dVmBlock2Djump_[2] = dVmBlock2Djump_[3] = ny*nz;
   dVmBlock2Djump_[4] = dVmBlock2Djump_[5] = ny*nz; //jump in x


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

   for(int ii=0;ii<27;ii++) offsetMap_[ii]=-100;
   offsetMap_[1*9+1*3+1]=ZZZ;
   offsetMap_[2*9+1*3+1]=PZZ;
   offsetMap_[1*9+2*3+1]=ZPZ;
   offsetMap_[0*9+1*3+1]=MZZ;
   offsetMap_[1*9+0*3+1]=ZMZ;
   offsetMap_[2*9+2*3+1]=PPZ;
   offsetMap_[0*9+2*3+1]=MPZ;
   offsetMap_[0*9+0*3+1]=MMZ;
   offsetMap_[2*9+0*3+1]=PMZ;
   offsetMap_[1*9+1*3+2]=ZZP;
   offsetMap_[1*9+1*3+0]=ZZM;
   offsetMap_[1*9+2*3+2]=ZPP;
   offsetMap_[1*9+2*3+0]=ZPM;
   offsetMap_[1*9+0*3+0]=ZMM;
   offsetMap_[1*9+0*3+2]=ZMP;
   offsetMap_[2*9+1*3+2]=PZP;
   offsetMap_[0*9+1*3+2]=MZP;
   offsetMap_[0*9+1*3+0]=MZM;
   offsetMap_[2*9+1*3+0]=PZM;

//   offsetsTest();
   precomputeCoefficients(anatomy);
   reorder_Coeff();

   
//   srand(1234);
//   double* VmM = VmBlock_.cBlock();
//   for(int ii=0;ii<nx*ny*nz;ii++)
//    { VmM[ii] = rand()/1000000.0;}
//
//   FGRDiff_simd_thread(1,nx-1,&(VmBlock_),dVmBlock_.cBlock());
//
//   for(int ii=1;ii<nx-1;ii++)
//   for(int jj=1;jj<ny-1;jj++)
//   for(int kk=1;kk<nz-1;kk++)
//   {
//     printf("(%d,%d,%d) : %f -> %f with %f\n",ii,jj,kk,VmBlock_(ii,jj,kk),dVmBlock_(ii,jj,kk), *((double*)&(diffCoefT2_(ii,jj,4*20*((int)(kk/4)) + 4* 5 + (kk%4)*2))) );
//   }
  
}

void FGRDiffusionOverlap::test()
{
  int tid = threadInfo_.teamRank();
  VectorDouble32 tVD32;
  //prepare Voltage
  double* BufferL = new double[nLocal_];
  double* BufferR = new double[nRemote_];
  unsigned nx = localGrid_.nx();
  unsigned ny = localGrid_.ny();
  unsigned nz = localGrid_.nz();

  for(int ii=0;ii<nLocal_;ii++) BufferL[ii]=ii;
  for(int ii=0;ii<nRemote_;ii++) BufferR[ii]=ii+1000000;

  printf("performing test\n");
  if(tid==0) reset_Coeff();

  //pass 1
  for(int ii=0;ii<nx*ny*nz;ii++) VmBlock_(ii)=0;
  for(int ii=0;ii<nx*ny*nz;ii++) dVmBlock_(ii)=0;
  L2_BarrierWithSync_Barrier(fgrBarrier_, &barrierHandle_[tid], threadInfo_.nThreads());
  updateLocalVoltage(BufferL);
  updateRemoteVoltageOld(BufferR);
  L2_BarrierWithSync_Barrier(fgrBarrier_, &barrierHandle_[tid], threadInfo_.nThreads());
  calc_overlap(tVD32);
  L2_BarrierWithSync_Barrier(fgrBarrier_, &barrierHandle_[tid], threadInfo_.nThreads());

  //pass 2
  Array3d<double> VmTmp(dVmBlock_);
  L2_BarrierWithSync_Barrier(fgrBarrier_, &barrierHandle_[tid], threadInfo_.nThreads());
  for(int ii=0;ii<nx*ny*nz;ii++) VmBlock_(ii)=0;
  for(int ii=0;ii<nx*ny*nz;ii++) dVmBlock_(ii)=0;
  L2_BarrierWithSync_Barrier(fgrBarrier_, &barrierHandle_[tid], threadInfo_.nThreads());
  updateLocalVoltage(BufferL);
  L2_BarrierWithSync_Barrier(fgrBarrier_, &barrierHandle_[tid], threadInfo_.nThreads());
  calc_overlap(tVD32);
  L2_BarrierWithSync_Barrier(fgrBarrier_, &barrierHandle_[tid], threadInfo_.nThreads());
  updateRemoteVoltage(BufferR);
  L2_BarrierWithSync_Barrier(fgrBarrier_, &barrierHandle_[tid], threadInfo_.nThreads());
  calc(tVD32);
  L2_BarrierWithSync_Barrier(fgrBarrier_, &barrierHandle_[tid], threadInfo_.nThreads());
  if(tid==0) compareVoltage(dVmBlock_,VmTmp);
  printf("end of test\n");
}

void FGRDiffusionOverlap::updateLocalVoltage(const double* VmLocal)
{
   startTimer(FGR_ArrayLocal2MatrixTimer);
   int tid = reactionThreadInfo_.teamRank();
   unsigned begin = localCopyOffset_[tid];
   unsigned end   = localCopyOffset_[tid+1];
   for (unsigned ii=begin; ii<end; ++ii)
   {
      int index = blockIndex_[ii];
      VmBlock_(index) = VmLocal[ii];
   }
   stopTimer(FGR_ArrayLocal2MatrixTimer);
}

void FGRDiffusionOverlap::updateRemoteVoltage(const double* VmRemote)
{
   startTimer(FGR_ArrayRemote2MatrixTimer);
   unsigned nx = localGrid_.nx();
   unsigned ny = localGrid_.ny();
   int tid = threadInfo_.teamRank();
   unsigned begin = remoteCopyOffset_[tid];
   unsigned end   = remoteCopyOffset_[tid+1];
   unsigned* bb = &blockIndex_[nLocal_];
   for (unsigned ii=begin; ii<end; ++ii)
   {
      int x = localTuple_[ii+nLocal_].x();
      int y = localTuple_[ii+nLocal_].y();
      int z = localTuple_[ii+nLocal_].z();

      if ( x == 0)         VmSlab_[0](0,y,z)=VmRemote[ii];
      else if( x == nx-1 ) VmSlab_[1](0,y,z)=VmRemote[ii];
      else if( y == 0 )    VmSlab_[2](0,x,z)=VmRemote[ii];
      else if( y == ny-1 ) VmSlab_[3](0,x,z)=VmRemote[ii];
      else if( z == BBzb ) VmSlab_[4](0,x,y)=VmRemote[ii];
      else if( z == BBze ) VmSlab_[5](0,x,y)=VmRemote[ii];
      else assert(0);
   }
   stopTimer(FGR_ArrayRemote2MatrixTimer);

}

void FGRDiffusionOverlap::updateRemoteVoltageOld(const double* VmRemote)
{
   startTimer(FGR_ArrayRemote2MatrixTimer);
   int tid = threadInfo_.teamRank();
   unsigned begin = remoteCopyOffset_[tid];
   unsigned end   = remoteCopyOffset_[tid+1];
   unsigned* bb = &blockIndex_[nLocal_];
   for (unsigned ii=begin; ii<end; ++ii)
   {
      int index = bb[ii];
      VmBlock_(index) = VmRemote[ii];
   }
   stopTimer(FGR_ArrayRemote2MatrixTimer);
}
//stencil on boundary
void FGRDiffusionOverlap::calc(VectorDouble32& dVm)
{
   int tid = threadInfo_.teamRank();

   #ifdef DIFF_TEST
   setAllVoltage(VmSlab_[0]);
   setAllVoltage(VmSlab_[1]);
   setAllVoltage(VmSlab_[2],1);
   setAllVoltage(VmSlab_[3],1);
   setAllVoltage(VmSlab_[4],2);
   setAllVoltage(VmSlab_[5],2);
   if(tid==0) printAllVoltage(VmSlab_[4],0);
   #endif

   startTimer(FGR_2D_StencilTimer);
   //how to divide the work among threads

   if (threadOffsetSimd_[tid] < threadOffsetSimd_[tid+1] )
   {
      FGRDiff_2D_z(0,threadOffset2D_[0][tid] ,threadOffset2D_[0][tid+1],dVmBlock_.cBlock());
      FGRDiff_2D_z(1,threadOffset2D_[0][tid] ,threadOffset2D_[0][tid+1],dVmBlock_.cBlock());
   }
   L2_BarrierWithSync_Barrier(fgrBarrier_, &barrierHandle_[tid], threadInfo_.nThreads());
   if (threadOffsetSimd_[tid] < threadOffsetSimd_[tid+1] )
   {
      FGRDiff_2D_z(2,threadOffset2D_[1][tid] ,threadOffset2D_[1][tid+1],dVmBlock_.cBlock());
      FGRDiff_2D_z(3,threadOffset2D_[1][tid] ,threadOffset2D_[1][tid+1],dVmBlock_.cBlock());
   }
   L2_BarrierWithSync_Barrier(fgrBarrier_, &barrierHandle_[tid], threadInfo_.nThreads());
   if (threadOffsetSimd_[tid] < threadOffsetSimd_[tid+1] )
   {
      FGRDiff_2D_xy(4,threadOffset2D_[2][tid] ,threadOffset2D_[2][tid+1]);
      FGRDiff_2D_xy(5,threadOffset2D_[2][tid] ,threadOffset2D_[2][tid+1]);
   }

   

   stopTimer(FGR_2D_StencilTimer);

   startTimer(FGR_Barrier2Timer);
   L2_BarrierWithSync_Barrier(fgrBarrier_, &barrierHandle_[tid], threadInfo_.nThreads());
   stopTimer(FGR_Barrier2Timer);

   #ifdef DIFF_TEST
   if(tid==0) printAllVoltage(dVmBlock_,2);
   Array3d<double> tmp(dVmBlock_);
   setAllVoltage(VmBlock_);
   copySlabToBlock();
   if (threadOffsetSimd_[tid] < threadOffsetSimd_[tid+1] )
      FGRDiff_simd_thread(threadOffsetSimd_[tid] ,threadOffsetSimd_[tid+1],&(VmBlock_),dVmBlock_.cBlock());
   L2_BarrierWithSync_Barrier(fgrBarrier_, &barrierHandle_[tid], threadInfo_.nThreads());
   if(tid==0) printAllVoltage(dVmBlock_,2);
   compareVoltage(dVmBlock_,tmp);
   #endif

   startTimer(FGR_Boundary2MatrixTimer);

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
   stopTimer(FGR_Boundary2MatrixTimer);
}

/** threaded simd version */
void FGRDiffusionOverlap::calc_overlap(VectorDouble32& dVm)
{
   int tid = threadInfo_.teamRank();

   Array3d<double> *VmTmp = &(VmBlock_);

   assert(VmBlock_.nz()%4 == 0); // make sure z is multiple of 4
   #ifdef DIFF_TEST
   setAllVoltage(*VmTmp);
   clearBlockBD();
   #endif
   
   startTimer(FGR_StencilTimer);
   if (threadOffsetSimd_[tid] < threadOffsetSimd_[tid+1] )
      FGRDiff_simd_thread(threadOffsetSimd_[tid] ,threadOffsetSimd_[tid+1],VmTmp,dVmBlock_.cBlock());
   stopTimer(FGR_StencilTimer);

   startTimer(FGR_Barrier2Timer);
   //L2_BarrierWithSync_Barrier(fgrBarrier_, &barrierHandle_[tid], threadInfo_.nThreads());
   stopTimer(FGR_Barrier2Timer);

}

//We're building the localTuple array.
void FGRDiffusionOverlap::buildTupleArray(const Anatomy& anatomy)
{
   localTuple_.resize(anatomy.size(), Tuple(0,0,0));
   BBzb=-1; //hope that it roles properly
   BBze=0;
   for (unsigned ii=0; ii<anatomy.size(); ++ii)
   {
      Tuple globalTuple = anatomy.globalTuple(ii);
      localTuple_[ii] = localGrid_.localTuple(globalTuple);
      if(localTuple_[ii].z()>BBze) BBze=localTuple_[ii].z();
      if(localTuple_[ii].z()<BBzb) BBzb=localTuple_[ii].z();
   }
   printf("BBzb=%d BBze=%d\n",BBzb,BBze);
}

void FGRDiffusionOverlap::buildBlockIndex(const Anatomy& anatomy)
{
   blockIndex_.resize(anatomy.size());
   for (unsigned ii=0; ii<anatomy.size(); ++ii)
   {
      Tuple globalTuple = anatomy.globalTuple(ii);
      Tuple ll = localGrid_.localTuple(globalTuple);
      blockIndex_[ii] = VmBlock_.tupleToIndex(ll.x(), ll.y(), ll.z());
   }
}

//void FGRDiffusionOverlap::precomputeCoefficientsPositive(const Anatomy& anatomy)
//{
//   unsigned nx = localGrid_.nx();
//   unsigned ny = localGrid_.ny();
//   unsigned nz = localGrid_.nz();
//   unsigned nxGlobal = anatomy.nx();
//   unsigned nyGlobal = anatomy.ny();
//   unsigned nzGlobal = anatomy.nz();
//   Vector hInv(1.0/anatomy.dx(), 1.0/anatomy.dy(), 1.0/anatomy.dz());
//   Vector h(anatomy.dx(), anatomy.dy(), anatomy.dz());
//   double gridCellVolume = h[0]*h[1]*h[2];
//   
//   SymmetricTensor sigmaZero = {0};
//   Array3d<SymmetricTensor> sigmaBlk(nx, ny, nz, sigmaZero);
//   Array3d<int> tissueBlk(nx, ny, nz, 0);
//
//   const vector<AnatomyCell>& cell = anatomy.cellArray();
//   for (unsigned ii=0; ii<anatomy.size(); ++ii)
//   {
//      unsigned ib = blockIndex_[ii];
//      sigmaBlk(ib) = anatomy.conductivity(ii);
//      tissueBlk(ib) = isTissue(anatomy.cellType(ii));
//   }
//
//   for (unsigned ii=0; ii<weight_.size(); ++ii)
//      for (unsigned jj=0; jj<19; ++jj)
//         weight_(ii).A[jj] = 0.0;
//
//   for (unsigned iCell=0; iCell<anatomy.nLocal(); ++iCell)
//   {
//      unsigned ib = blockIndex_[iCell];
//      int have2D=0;
//      //check if any quadrant is formed.
//      for(int dir=0;dir<3;dir++) //z,x,y
//      {
//        //(110) (011) (101)
//        for(int ii=-1;ii<2;ii+=2)
//        for(int jj=-1;jj<2;jj+=2)
//        {
//          int have4=0;
//          int coordi[3]={0,0,0};
//          SymmetricTensor avgSigma ={0};
//          for(int ll=0;ll<2;ll++)
//          for(int mm=0;mm<2;ll++)
//          {
//            int i1=dir;
//            int i2=(dir+1)%3;
//            coordi[i1]=ll*ii;
//            coordi[i2]=mm*jj;
//
//            offset = sigmaBlk.tupleToIndex(coordi[0]+1,coordi[1]+1,coordi[2]+1) - sigmaBlk.tupleToIndex(1,1,1);
//            if (tissueBlk(ib + offset) == 0) break;
//
//            avgSigma = avgSigma + sigmaBlk(ib+offset);
//            have4++; 
//          }
//          if ( have4 == 4 ) 
//          {
//            int coordi[3]={0,0,0};
//            int signA[3]={0,0,0};
//            for(int ll=0;ll<2;ll++)
//            for(int mm=0;mm<2;ll++)
//            {
//              int i1=dir;
//              int i2=(dir+1)%3;
//              coordi[i1]=ll*ii;
//              coordi[i2]=mm*jj;
//
//              signA[i1]=ii*(2*ll-1);
//              signA[i2]=jj*(2*mm-1);
// 
//              Vector signV(signA[0],signA[1],signA[2]);
//              Vector weightV = avgSigma * signV;
//
//              offset = 9*(coordi[0]+1)+3*(coordi[1]+1)+(coordi[2]+1);
//              assert(offsetMap_[offset] != -100);
//              weight_(ib).A[offsetMap_[offset]] += weightV.sum();
//            }
//            have2D++;
//          }
//        }         
//      }
//      if ( have2D == 0)
//      {
//         //check if 1D 
//      }
//      
//      double sum = 0;
//      for (unsigned ii=0; ii<19; ++ii)
//         sum += weight_(ib).A[ii];
//      assert(abs(sum) < weightSumTolerance);
//   }
////   printAllWeights(tissueBlk);
//}


void FGRDiffusionOverlap::precomputeCoefficients(const Anatomy& anatomy)
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

      double sum = weight_(ib).A[0];
      for (unsigned ii=1; ii<19; ++ii)
      {
         sum += weight_(ib).A[ii];
         A0_(ib) -= weight_(ib).A[ii];
      }
      assert(abs(sum) < weightSumTolerance);
   }
   //printAllWeights(tissueBlk);
}

void FGRDiffusionOverlap::mkTissueArray(
   const Array3d<int>& tissueBlk, int ib, int* tissue)
{
   for (unsigned ii=0; ii<19; ++ii)
      tissue[ii] = tissueBlk(ib + offset_[ii]);
}


Vector FGRDiffusionOverlap::f1(int ib, int iFace, const Vector& h,
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


void FGRDiffusionOverlap::printAllWeights(const Array3d<int>& tissue)
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

void FGRDiffusionOverlap::printAllVoltage(Array3d<double>& Voltage,int map)
{
   int ll[3]={Voltage.nx(),Voltage.ny(),Voltage.nz()};
   int xi[3];
   for(xi[0]=0;xi[0]<ll[map];xi[0]++)
   for(xi[1]=0;xi[1]<ll[(map+1)%3];xi[1]++)
   for(xi[2]=0;xi[2]<ll[(map+2)%3];xi[2]++)
   {
      int ei[3]={xi[(map*2)%3],xi[(map*2+1)%3],xi[(map*2+2)%3]};
      printf("Voltage: %5d %5d %5d %20.12e\n",ei[0],ei[1],ei[2],Voltage(ei[0],ei[1],ei[2]));
   }
}

void FGRDiffusionOverlap::setAllVoltage(Array3d<double>& Voltage, int cut)
{
   for(int xx=0;xx<Voltage.nx();xx++)
   for(int yy=0;yy<Voltage.ny();yy++)
   for(int zz=0;zz<Voltage.nz();zz++)
   {
     Voltage(xx,yy,zz)=zz + 100* ( yy + 100* xx );
     if (cut>0)
     {
       if( yy==0) Voltage(xx,yy,zz)=0;
       if( yy==Voltage.ny()-1) Voltage(xx,yy,zz)=0;
     }
     if (cut>1)
     {
       if( zz==0) Voltage(xx,yy,zz)=0;
       if( zz==Voltage.nz()-1) Voltage(xx,yy,zz)=0;
     }
   }
}

void FGRDiffusionOverlap::reorder_Coeff()
{
  uint32_t idx=0,ll,ii,jj,kk;
  const uint32_t Nx2 = localGrid_.nx();
  const uint32_t Ny2 = localGrid_.ny();
  const uint32_t Nz2 = localGrid_.nz();

  diffCoefT2_.resize(Nx2,Ny2,Nz2*20,0.0);
  diffCoefT3_[0].resize(1,Ny2,Nz2*5,0.0);  //yz
  diffCoefT3_[1].resize(1,Ny2,Nz2*5,0.0);  

  diffCoefT3_[2].resize(1,Nx2,Nz2*5,0.0);  //xz
  diffCoefT3_[3].resize(1,Nx2,Nz2*5,0.0);

  diffCoefT3_[4].resize(1,Nx2,((int)(Ny2+3)/4)*4*5,0.0);  //xy
  diffCoefT3_[5].resize(1,Nx2,((int)(Ny2+3)/4)*4*5,0.0);

  for(int ii=0;ii<nLocal_;ii++)
  {
    int xx = localTuple_[ii].x();
    int yy = localTuple_[ii].y();
    int zz = localTuple_[ii].z();
    int xx4 = (int)(xx/4);
    int yy4 = (int)(yy/4);
    int zz4 = (int)(zz/4);
    int xp = xx + 1;
    int yp = yy + 1;
    int zp = zz + 1;
    int xm = xx - 1;
    int ym = yy - 1;
    int zm = zz - 1;
    int xp4 = (int)(xp/4);
    int yp4 = (int)(yp/4);
    int zp4 = (int)(zp/4);
    int xm4 = (int)(xm/4);
    int ym4 = (int)(ym/4);
    int zm4 = (int)(zm/4);

    assert(zz > 0);
    assert(zz < Nz2-1);

    diffCoefT2_(xx,yy,4*20*zp4 + 4* 0 + zp%4 ) =weight_(xx,yy,zz).A[ZZP];
    diffCoefT2_(xx,yy,4*20*zp4 + 4* 1 + zp%4 ) =weight_(xx,yy,zz).A[ZMP];
    diffCoefT2_(xx,yy,4*20*zp4 + 4* 2 + zp%4 ) =weight_(xx,yy,zz).A[MZP];
    diffCoefT2_(xx,yy,4*20*zp4 + 4* 3 + zp%4 ) =weight_(xx,yy,zz).A[ZPP];
    diffCoefT2_(xx,yy,4*20*zp4 + 4* 6 + zp%4 ) =weight_(xx,yy,zz).A[PZP];
                                 
    *((double*)&(diffCoefT2_(xx,yy,4*20*zz4 + 4* 4 + (zz%4)*2))) = A0_(xx,yy,zz);
//    diffCoefT2_(xx,yy,4*20*zz4 + 4* 5 + zz%4 ) =weight_(xx,yy,zz).A[ZZZ]=0;
    diffCoefT2_(xx,yy,4*20*zz4 + 4* 7 + zz%4 ) =weight_(xx,yy,zz).A[PMZ];
    diffCoefT2_(xx,yy,4*20*zz4 + 4* 8 + zz%4 ) =weight_(xx,yy,zz).A[MMZ];
    diffCoefT2_(xx,yy,4*20*zz4 + 4* 9 + zz%4 ) =weight_(xx,yy,zz).A[MPZ];
    diffCoefT2_(xx,yy,4*20*zz4 + 4*10 + zz%4 ) =weight_(xx,yy,zz).A[PPZ];
    diffCoefT2_(xx,yy,4*20*zz4 + 4*11 + zz%4 ) =weight_(xx,yy,zz).A[ZMZ];
    diffCoefT2_(xx,yy,4*20*zz4 + 4*12 + zz%4 ) =weight_(xx,yy,zz).A[MZZ];
    diffCoefT2_(xx,yy,4*20*zz4 + 4*13 + zz%4 ) =weight_(xx,yy,zz).A[ZPZ];
    diffCoefT2_(xx,yy,4*20*zz4 + 4*14 + zz%4 ) =weight_(xx,yy,zz).A[PZZ];

    diffCoefT2_(xx,yy,4*20*zm4 + 4*15 + zm%4 ) =weight_(xx,yy,zz).A[ZZM];
    diffCoefT2_(xx,yy,4*20*zm4 + 4*16 + zm%4 ) =weight_(xx,yy,zz).A[ZMM];
    diffCoefT2_(xx,yy,4*20*zm4 + 4*17 + zm%4 ) =weight_(xx,yy,zz).A[MZM];
    diffCoefT2_(xx,yy,4*20*zm4 + 4*18 + zm%4 ) =weight_(xx,yy,zz).A[ZPM];
    diffCoefT2_(xx,yy,4*20*zm4 + 4*19 + zm%4 ) =weight_(xx,yy,zz).A[PZM];

    if (xx == 1)
    {
      diffCoefT3_[0](0,yy,4*5*zz4 + 4*0 + zz%4) = weight_(xx,yy,zz).A[MZZ];
      diffCoefT3_[0](0,yy,4*5*zz4 + 4*1 + zz%4) = weight_(xx,yy,zz).A[MPZ];
      diffCoefT3_[0](0,yy,4*5*zz4 + 4*2 + zz%4) = weight_(xx,yy,zz).A[MMZ];
      diffCoefT3_[0](0,yy,4*5*zp4 + 4*3 + zp%4) = weight_(xx,yy,zz).A[MZP];
      diffCoefT3_[0](0,yy,4*5*zm4 + 4*4 + zm%4) = weight_(xx,yy,zz).A[MZM];
    }
    if (xx == Nx2-2)
    {
      diffCoefT3_[1](0,yy,4*5*zz4 + 4*0 + zz%4) = weight_(xx,yy,zz).A[PZZ];
      diffCoefT3_[1](0,yy,4*5*zz4 + 4*1 + zz%4) = weight_(xx,yy,zz).A[PPZ];
      diffCoefT3_[1](0,yy,4*5*zz4 + 4*2 + zz%4) = weight_(xx,yy,zz).A[PMZ];
      diffCoefT3_[1](0,yy,4*5*zp4 + 4*3 + zp%4) = weight_(xx,yy,zz).A[PZP];
      diffCoefT3_[1](0,yy,4*5*zm4 + 4*4 + zm%4) = weight_(xx,yy,zz).A[PZM];
    }
    if( yy == 1 )
    {
      diffCoefT3_[2](0,xx,4*5*zz4 + 4*0 + zz%4) = weight_(xx,yy,zz).A[ZMZ];
      diffCoefT3_[2](0,xx,4*5*zz4 + 4*1 + zz%4) = weight_(xx,yy,zz).A[PMZ];
      diffCoefT3_[2](0,xx,4*5*zz4 + 4*2 + zz%4) = weight_(xx,yy,zz).A[MMZ];
      diffCoefT3_[2](0,xx,4*5*zp4 + 4*3 + zp%4) = weight_(xx,yy,zz).A[ZMP];
      diffCoefT3_[2](0,xx,4*5*zm4 + 4*4 + zm%4) = weight_(xx,yy,zz).A[ZMM];
    }
    if( yy == Ny2-2 )
    {
      diffCoefT3_[3](0,xx,4*5*zz4 + 4*0 + zz%4) = weight_(xx,yy,zz).A[ZPZ];
      diffCoefT3_[3](0,xx,4*5*zz4 + 4*1 + zz%4) = weight_(xx,yy,zz).A[PPZ];
      diffCoefT3_[3](0,xx,4*5*zz4 + 4*2 + zz%4) = weight_(xx,yy,zz).A[MPZ];
      diffCoefT3_[3](0,xx,4*5*zp4 + 4*3 + zp%4) = weight_(xx,yy,zz).A[ZPP];
      diffCoefT3_[3](0,xx,4*5*zm4 + 4*4 + zm%4) = weight_(xx,yy,zz).A[ZPM];
    }
    if( zz == BBzb+1 )
    {
      diffCoefT3_[4](0,xx,4*5*yy4 + 4*0 + yy%4) = weight_(xx,yy,zz).A[ZZM];
      diffCoefT3_[4](0,xx,4*5*yy4 + 4*1 + yy%4) = weight_(xx,yy,zz).A[PZM];
      diffCoefT3_[4](0,xx,4*5*yy4 + 4*2 + yy%4) = weight_(xx,yy,zz).A[MZM];
      diffCoefT3_[4](0,xx,4*5*yp4 + 4*3 + yp%4) = weight_(xx,yy,zz).A[ZPM];
      diffCoefT3_[4](0,xx,4*5*ym4 + 4*4 + ym%4) = weight_(xx,yy,zz).A[ZMM];
    }
    if( zz == BBze-1 )
    {
      diffCoefT3_[5](0,xx,4*5*yy4 + 4*0 + yy%4) = weight_(xx,yy,zz).A[ZZP];
      diffCoefT3_[5](0,xx,4*5*yy4 + 4*1 + yy%4) = weight_(xx,yy,zz).A[PZP];
      diffCoefT3_[5](0,xx,4*5*yy4 + 4*2 + yy%4) = weight_(xx,yy,zz).A[MZP];
      diffCoefT3_[5](0,xx,4*5*yp4 + 4*3 + yp%4) = weight_(xx,yy,zz).A[ZPP];
      diffCoefT3_[5](0,xx,4*5*ym4 + 4*4 + ym%4) = weight_(xx,yy,zz).A[ZMP];
    }
  }
}

void FGRDiffusionOverlap::reset_Coeff()
{
  uint32_t idx=0,ll,ii,jj,kk;
  const uint32_t Nx2 = localGrid_.nx();
  const uint32_t Ny2 = localGrid_.ny();
  const uint32_t Nz2 = localGrid_.nz();

  for(int ii=0;ii<nLocal_;ii++)
  {
    int xx = localTuple_[ii].x();
    int yy = localTuple_[ii].y();
    int zz = localTuple_[ii].z();
    int xx4 = (int)(xx/4);
    int yy4 = (int)(yy/4);
    int zz4 = (int)(zz/4);
    int xp = xx + 1;
    int yp = yy + 1;
    int zp = zz + 1;
    int xm = xx - 1;
    int ym = yy - 1;
    int zm = zz - 1;
    int xp4 = (int)(xp/4);
    int yp4 = (int)(yp/4);
    int zp4 = (int)(zp/4);
    int xm4 = (int)(xm/4);
    int ym4 = (int)(ym/4);
    int zm4 = (int)(zm/4);

    assert(zz > 0);
    assert(zz < Nz2-1);

    diffCoefT2_(xx,yy,4*20*zp4 + 4* 0 + zp%4 ) =weight_(xx,yy,zz).A[ZZP]=1;
    diffCoefT2_(xx,yy,4*20*zp4 + 4* 1 + zp%4 ) =weight_(xx,yy,zz).A[ZMP]=1;
    diffCoefT2_(xx,yy,4*20*zp4 + 4* 2 + zp%4 ) =weight_(xx,yy,zz).A[MZP]=1;
    diffCoefT2_(xx,yy,4*20*zp4 + 4* 3 + zp%4 ) =weight_(xx,yy,zz).A[ZPP]=1;
    diffCoefT2_(xx,yy,4*20*zp4 + 4* 6 + zp%4 ) =weight_(xx,yy,zz).A[PZP]=1;
                                 
    *((double*)&(diffCoefT2_(xx,yy,4*20*zz4 + 4* 4 + (zz%4)*2))) = A0_(xx,yy,zz)=1;
//    diffCoefT2_(xx,yy,4*20*zz4 + 4* 5 + zz%4 ) =weight_(xx,yy,zz).A[ZZZ]=0=1;
    diffCoefT2_(xx,yy,4*20*zz4 + 4* 7 + zz%4 ) =weight_(xx,yy,zz).A[PMZ]=1;
    diffCoefT2_(xx,yy,4*20*zz4 + 4* 8 + zz%4 ) =weight_(xx,yy,zz).A[MMZ]=1;
    diffCoefT2_(xx,yy,4*20*zz4 + 4* 9 + zz%4 ) =weight_(xx,yy,zz).A[MPZ]=1;
    diffCoefT2_(xx,yy,4*20*zz4 + 4*10 + zz%4 ) =weight_(xx,yy,zz).A[PPZ]=1;
    diffCoefT2_(xx,yy,4*20*zz4 + 4*11 + zz%4 ) =weight_(xx,yy,zz).A[ZMZ]=1;
    diffCoefT2_(xx,yy,4*20*zz4 + 4*12 + zz%4 ) =weight_(xx,yy,zz).A[MZZ]=1;
    diffCoefT2_(xx,yy,4*20*zz4 + 4*13 + zz%4 ) =weight_(xx,yy,zz).A[ZPZ]=1;
    diffCoefT2_(xx,yy,4*20*zz4 + 4*14 + zz%4 ) =weight_(xx,yy,zz).A[PZZ]=1;

    diffCoefT2_(xx,yy,4*20*zm4 + 4*15 + zm%4 ) =weight_(xx,yy,zz).A[ZZM]=1;
    diffCoefT2_(xx,yy,4*20*zm4 + 4*16 + zm%4 ) =weight_(xx,yy,zz).A[ZMM]=1;
    diffCoefT2_(xx,yy,4*20*zm4 + 4*17 + zm%4 ) =weight_(xx,yy,zz).A[MZM]=1;
    diffCoefT2_(xx,yy,4*20*zm4 + 4*18 + zm%4 ) =weight_(xx,yy,zz).A[ZPM]=1;
    diffCoefT2_(xx,yy,4*20*zm4 + 4*19 + zm%4 ) =weight_(xx,yy,zz).A[PZM]=1;

    if (xx == 1)
    {
      diffCoefT3_[0](0,yy,4*5*zz4 + 4*0 + zz%4) = weight_(xx,yy,zz).A[MZZ];
      diffCoefT3_[0](0,yy,4*5*zz4 + 4*1 + zz%4) = weight_(xx,yy,zz).A[MPZ];
      diffCoefT3_[0](0,yy,4*5*zz4 + 4*2 + zz%4) = weight_(xx,yy,zz).A[MMZ];
      diffCoefT3_[0](0,yy,4*5*zp4 + 4*3 + zp%4) = weight_(xx,yy,zz).A[MZP];
      diffCoefT3_[0](0,yy,4*5*zm4 + 4*4 + zm%4) = weight_(xx,yy,zz).A[MZM];
    }
    if (xx == Nx2-2)
    {
      diffCoefT3_[1](0,yy,4*5*zz4 + 4*0 + zz%4) = weight_(xx,yy,zz).A[PZZ];
      diffCoefT3_[1](0,yy,4*5*zz4 + 4*1 + zz%4) = weight_(xx,yy,zz).A[PPZ];
      diffCoefT3_[1](0,yy,4*5*zz4 + 4*2 + zz%4) = weight_(xx,yy,zz).A[PMZ];
      diffCoefT3_[1](0,yy,4*5*zp4 + 4*3 + zp%4) = weight_(xx,yy,zz).A[PZP];
      diffCoefT3_[1](0,yy,4*5*zm4 + 4*4 + zm%4) = weight_(xx,yy,zz).A[PZM];
    }
    if( yy == 1 )
    {
      diffCoefT3_[2](0,xx,4*5*zz4 + 4*0 + zz%4) = weight_(xx,yy,zz).A[ZMZ];
      diffCoefT3_[2](0,xx,4*5*zz4 + 4*1 + zz%4) = weight_(xx,yy,zz).A[PMZ];
      diffCoefT3_[2](0,xx,4*5*zz4 + 4*2 + zz%4) = weight_(xx,yy,zz).A[MMZ];
      diffCoefT3_[2](0,xx,4*5*zp4 + 4*3 + zp%4) = weight_(xx,yy,zz).A[ZMP];
      diffCoefT3_[2](0,xx,4*5*zm4 + 4*4 + zm%4) = weight_(xx,yy,zz).A[ZMM];
    }
    if( yy == Ny2-2 )
    {
      diffCoefT3_[3](0,xx,4*5*zz4 + 4*0 + zz%4) = weight_(xx,yy,zz).A[ZPZ];
      diffCoefT3_[3](0,xx,4*5*zz4 + 4*1 + zz%4) = weight_(xx,yy,zz).A[PPZ];
      diffCoefT3_[3](0,xx,4*5*zz4 + 4*2 + zz%4) = weight_(xx,yy,zz).A[MPZ];
      diffCoefT3_[3](0,xx,4*5*zp4 + 4*3 + zp%4) = weight_(xx,yy,zz).A[ZPP];
      diffCoefT3_[3](0,xx,4*5*zm4 + 4*4 + zm%4) = weight_(xx,yy,zz).A[ZPM];
    }
    if( zz == BBzb+1 )
    {
      diffCoefT3_[4](0,xx,4*5*yy4 + 4*0 + yy%4) = weight_(xx,yy,zz).A[ZZM];
      diffCoefT3_[4](0,xx,4*5*yy4 + 4*1 + yy%4) = weight_(xx,yy,zz).A[PZM];
      diffCoefT3_[4](0,xx,4*5*yy4 + 4*2 + yy%4) = weight_(xx,yy,zz).A[MZM];
      diffCoefT3_[4](0,xx,4*5*yp4 + 4*3 + yp%4) = weight_(xx,yy,zz).A[ZPM];
      diffCoefT3_[4](0,xx,4*5*ym4 + 4*4 + ym%4) = weight_(xx,yy,zz).A[ZMM];
    }
    if( zz == BBze-1 )
    {
      diffCoefT3_[5](0,xx,4*5*yy4 + 4*0 + yy%4) = weight_(xx,yy,zz).A[ZZP];
      diffCoefT3_[5](0,xx,4*5*yy4 + 4*1 + yy%4) = weight_(xx,yy,zz).A[PZP];
      diffCoefT3_[5](0,xx,4*5*yy4 + 4*2 + yy%4) = weight_(xx,yy,zz).A[MZP];
      diffCoefT3_[5](0,xx,4*5*yp4 + 4*3 + yp%4) = weight_(xx,yy,zz).A[ZPP];
      diffCoefT3_[5](0,xx,4*5*ym4 + 4*4 + ym%4) = weight_(xx,yy,zz).A[ZMP];
    }
  }
}
// simdiazed version
// 'start' cannot be in the middle. 'start' must point to (n,m,0). 
void
FGRDiffusionOverlap::FGRDiff_simd_thread(const uint32_t bx,const int32_t ex, Array3d<double>* VmTmp, double* out0)
{
  int ii;

  const unsigned Nx2 = VmTmp->nx();
  const unsigned Ny2 = VmTmp->ny();
  const unsigned Nz2 = VmTmp->nz();

  double* VmM = VmTmp->cBlock();
  double* out;
  WeightType *diffC = diffCoefT2_.cBlock();

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
  vector4double double_vec;
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
 
    WeightType *simd_diff_ = diffC + start*20; 
 
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
        x =  vec_madd(vec_ld(4*6*WTSZ,simd_diff_)  , my_xp1_y_z, \
             vec_madd(vec_ld(4*3*WTSZ,simd_diff_)  , my_x_yp1_z, \
             vec_madd(vec_ld(4*2*WTSZ,simd_diff_)  , my_xm1_y_z, \
             vec_madd(vec_ld(4*1*WTSZ,simd_diff_)  , my_x_ym1_z, \
             vec_mul( vec_ld(4*0*WTSZ,simd_diff_)  , my_x_y_z)))));
 
    #define calc_zz(x) \
        double_vec = vec_ld(4*4*WTSZ,(double*)simd_diff_); \
        x = vec_madd( vec_ld(4*14*WTSZ,simd_diff_)  , my_xp1_y_z, \
            vec_madd( vec_ld(4*13*WTSZ,simd_diff_)  , my_x_yp1_z, \
            vec_madd( vec_ld(4*12*WTSZ,simd_diff_)  , my_xm1_y_z, \
            vec_madd( vec_ld(4*11*WTSZ,simd_diff_)  , my_x_ym1_z, \
            vec_madd( vec_ld(4*10*WTSZ,simd_diff_)  , my_xp1_yp1_z, \
            vec_madd( vec_ld(4*9*WTSZ,simd_diff_)  , my_xm1_yp1_z, \
            vec_madd( vec_ld(4*8*WTSZ,simd_diff_)  , my_xm1_ym1_z, \
            vec_madd( vec_ld(4*7*WTSZ,simd_diff_)  , my_xp1_ym1_z, \
            vec_mul (double_vec  , my_x_y_z)))))))));
//            vec_mul ( vec_ld(4*5*WTSZ,(double*)simd_diff_)  , my_x_y_z)))))))));
            //vec_mul ( vec_ld(4*5*WTSZ,simd_diff_)  , my_x_y_z)))))))));
 
    #define calc_zm(x) \
        x = vec_madd( vec_ld(4*19*WTSZ,simd_diff_)  , my_xp1_y_z, \
            vec_madd( vec_ld(4*18*WTSZ,simd_diff_)  , my_x_yp1_z, \
            vec_madd( vec_ld(4*17*WTSZ,simd_diff_)  , my_xm1_y_z, \
            vec_madd( vec_ld(4*16*WTSZ,simd_diff_)  , my_x_ym1_z, \
            vec_mul ( vec_ld(4*15*WTSZ,simd_diff_)  , my_x_y_z)))));
 
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
 
    double* simd_diff2 = (double*)simd_diff_;
    load_my_vectors;
 
    calc_zp(B0);
    calc_zz(Sum1);
    calc_zm(B2);
 
    Sum2 = vec_sldw(my_zero_vec,B2,3);
 
//    cout << "coeff:" << *(simd_diff2 + 2*5 + 1) <<  " " <<  *(simd_diff2 + 2*5 + 2) <<  " " << *(simd_diff2 + 2*5 + 3) <<  " " << endl;

//    cout << "sum1:" << Sum1[0] << endl;

    for (ii=0;ii<chunk_size;ii+=4)
    {

      simd_diff_ += 20*4;



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

//      printf("Sum1:%f %f %f %f\n",Sum1[0],Sum1[1],Sum1[2],Sum1[3]);
//      printf("double_vec:%f %f %f %f\n",double_vec[0],double_vec[1],double_vec[2],double_vec[3]);
//      printf("double_vec2:%f %f %f %f\n",*((double*)(simd_diff_+4*5)),*((double*)(simd_diff_+4*5+2)),*((double*)(simd_diff_+4*5+4)),*((double*)(simd_diff_+4*5+6)));
    }
  }
};

void
FGRDiffusionOverlap::FGRDiff_2D_xy(uint32_t slabID, const uint32_t by,const int32_t ey)
{
#ifdef BGQ 
  int ii;

  const unsigned Ny2 = VmSlab_[slabID].ny(); //ny means the second dimension
  const unsigned Nz2 = VmSlab_[slabID].nz(); //nz means the third dimension
  const unsigned out_jump = dVmBlock2Djump_[slabID];
  const unsigned ib0 = 0;
  const unsigned ib1 = dVmBlock_.nz()*8;
  const unsigned ib2 = ib1*2;
  const unsigned ib3 = ib1*3;
  const unsigned ib4 = ib1/2;

  double* VmM = VmSlab_[slabID].cBlock();
  double* out;
  double* out0 = dVmBlock_.cBlock();
  out0+=dVmBlock2Doffset_[slabID] + by*out_jump;
  WeightType *diffC = diffCoefT3_[slabID].cBlock();

  int xym1z_ =   ((0 ) *Ny2 + (-1)) * Nz2;
  int xyp1z_ =   ((0 ) *Ny2 + (+1)) * Nz2;

  vector4double B0,Sum1,B2,C0,C1,C2,Sum0,Sum2,Sum;
  vector4double my_0        ;
  vector4double my_x_y_z    ;
  vector4double my_x_yp1_z  ;
  vector4double my_x_ym1_z  ;
  vector4double my_zero_vec = vec_splats(0.0);
//  vector4double my_one_vec = vec_splats(0.1);
  vector4double y1,y2,y3;

  for(int yy=by;yy<ey;yy++)
  {
    int start = VmSlab_[slabID].tupleToIndex(0,yy,0);
    int end   = VmSlab_[slabID].tupleToIndex(0,yy,Nz2);
    int chunk_size = end - start;

    double* phi_x_ym1_z   = VmM + start + xym1z_;
    double* phi_x_y_z     = VmM + start ;
    double* phi_x_yp1_z   = VmM + start + xyp1z_;
 
    out = out0;
 
    WeightType *simd_diff_ = diffC + start*5; 
 
    #define load_my_0 \
	asm("lfdx %0,%1,%2" : "=v"(my_0) : "b"(ib0) , "b"(out)); \
	asm("lfdx %0,%1,%2" : "=v"(y1) : "b"(ib1) , "b"(out)); \
	my_0 = vec_sldw(my_0,y1,1); \
	asm("lfdx %0,%1,%2" : "=v"(y2) : "b"(ib2) , "b"(out)); \
	my_0 = vec_sldw(my_0,y2,1); \
	asm("lfdx %0,%1,%2" : "=v"(y3) : "b"(ib3) , "b"(out)); \
	my_0 = vec_sldw(my_0,y3,1);

    #define scatter_store(x) \
      asm("stfdx %0,%1,%2" : : "v"(x) , "b"(ib0) , "b"(out)); \
      y1 = vec_sldw(x,x,1); \
      asm("stfdx %0,%1,%2" : : "v"(y1) , "b"(ib1) , "b"(out)); \
      y2 = vec_sldw(x,x,2); \
      asm("stfdx %0,%1,%2" : : "v"(y2) , "b"(ib2) , "b"(out)); \
      y3 = vec_sldw(x,x,3); \
      asm("stfdx %0,%1,%2" : : "v"(y3) , "b"(ib3) , "b"(out));

    #define load_my_vectors_2D  \
        my_x_y_z      =vec_ld(0,      phi_x_y_z   );\
        my_x_yp1_z    =vec_ld(0,      phi_x_yp1_z );\
        my_x_ym1_z    =vec_ld(0,      phi_x_ym1_z );
 
    #define calc_zp_2D(x) \
        x = vec_mul( vec_ld(4*3*WTSZ,simd_diff_)  , my_x_y_z);
 
    #define calc_zz_2D(x) \
        x = vec_madd( vec_ld(4*2*WTSZ,simd_diff_)  , my_x_ym1_z, \
            vec_madd( vec_ld(4*1*WTSZ,simd_diff_)  , my_x_yp1_z, \
            vec_madd( vec_ld(4*0*WTSZ,simd_diff_)  , my_x_y_z, my_0))); \
 
    #define calc_zm_2D(x) \
        x = vec_mul( vec_ld(4*4*WTSZ,simd_diff_), my_x_y_z);
 
    #define shift_pointers_2D \
            phi_x_yp1_z   +=4; \
            phi_x_ym1_z   +=4; \
            phi_x_y_z     +=4;
 
    load_my_vectors_2D;
    load_my_0;
 
    calc_zp_2D(B0);
    calc_zz_2D(Sum1);
    calc_zm_2D(B2);
 
    Sum2 = vec_sldw(my_zero_vec,B2,3);
 
//    cout << "sum1:" << Sum1[0] << endl;

    for (ii=0;ii<chunk_size;ii+=4)
    {

      simd_diff_ += 5*4;

      shift_pointers_2D; //note that these pointers one step further
      load_my_vectors_2D;
 
      calc_zp_2D(C0);
      Sum0 = vec_sldw(B0,C0,1);
      B0=C0;
 
      Sum = vec_add(vec_add(Sum2,Sum1),Sum0);
      scatter_store(Sum);
      out+=ib4;
 
      load_my_0;
      calc_zz_2D(Sum1);
      calc_zm_2D(C2);
      Sum2 = vec_sldw(B2,C2,3);
      B2 = C2;

//      printf("Sum1:%f %f %f %f\n",Sum1[0],Sum1[1],Sum1[2],Sum1[3]);
//      printf("double_vec:%f %f %f %f\n",double_vec[0],double_vec[1],double_vec[2],double_vec[3]);
//      printf("double_vec2:%f %f %f %f\n",*((double*)(simd_diff_+4*5)),*((double*)(simd_diff_+4*5+2)),*((double*)(simd_diff_+4*5+4)),*((double*)(simd_diff_+4*5+6)));
    }
    out0 += out_jump;
  }
  #undef load_my_vectors_2D  
  #undef calc_zp_2D
  #undef calc_zz_2D
  #undef calc_zm_2D 
  #undef shift_pointers_2D 
#endif
}

void
FGRDiffusionOverlap::FGRDiff_2D_z(uint32_t slabID, const uint32_t by,const int32_t ey, double* out0)
{
  int ii;

  const unsigned Ny2 = VmSlab_[slabID].ny(); //ny means the second dimension
  const unsigned Nz2 = VmSlab_[slabID].nz(); //nz means the third dimension
  const unsigned out_jump = dVmBlock2Djump_[slabID];

  double* VmM = VmSlab_[slabID].cBlock();
  double* out;
  out0 += dVmBlock2Doffset_[slabID] + by*dVmBlock2Djump_[slabID];
  WeightType *diffC = diffCoefT3_[slabID].cBlock();

  int xym1z_ =   ((0 ) *Ny2 + (-1)) * Nz2;
  int xyp1z_ =   ((0 ) *Ny2 + (+1)) * Nz2;

  vector4double B0,Sum1,B2,C0,C1,C2,Sum0,Sum2,Sum;
  vector4double my_0        ;
  vector4double my_x_y_z    ;
  vector4double my_x_yp1_z  ;
  vector4double my_x_ym1_z  ;
  vector4double my_zero_vec = vec_splats(0.0);

  for(int yy=by;yy<ey;yy++)
  {
    int start = VmSlab_[slabID].tupleToIndex(0,yy,0);
    int end   = VmSlab_[slabID].tupleToIndex(0,yy,Nz2);
    int chunk_size = end - start;

    double* phi_x_ym1_z   = VmM + start + xym1z_;
    double* phi_x_y_z     = VmM + start ;
    double* phi_x_yp1_z   = VmM + start + xyp1z_;
 
    out = out0;
 
    WeightType *simd_diff_ = diffC + start*5; 
 
    #define load_my_vectors_2D  \
        my_x_y_z      =vec_ld(0,      phi_x_y_z   );\
        my_x_yp1_z    =vec_ld(0,      phi_x_yp1_z );\
        my_x_ym1_z    =vec_ld(0,      phi_x_ym1_z );
 
    #define calc_zp_2D(x) \
        x = vec_mul( vec_ld(4*3*WTSZ,simd_diff_)  , my_x_y_z);
 
    #define calc_zz_2D(x) \
        x = vec_madd( vec_ld(4*2*WTSZ,simd_diff_)  , my_x_ym1_z, \
            vec_madd( vec_ld(4*1*WTSZ,simd_diff_)  , my_x_yp1_z, \
            vec_madd( vec_ld(4*0*WTSZ,simd_diff_)  , my_x_y_z, my_0))); \
 
    #define calc_zm_2D(x) \
        x = vec_mul( vec_ld(4*4*WTSZ,simd_diff_), my_x_y_z);
 
    #define shift_pointers_2D \
            phi_x_yp1_z   +=4; \
            phi_x_ym1_z   +=4; \
            phi_x_y_z     +=4;
 
    double* simd_diff2 = (double*)simd_diff_;
    load_my_vectors_2D;
    my_0=vec_ld(0,out);
 
    calc_zp_2D(B0);
    calc_zz_2D(Sum1);
    calc_zm_2D(B2);
 
    Sum2 = vec_sldw(my_zero_vec,B2,3);
 
//    cout << "coeff:" << *(simd_diff2 + 2*5 + 1) <<  " " <<  *(simd_diff2 + 2*5 + 2) <<  " " << *(simd_diff2 + 2*5 + 3) <<  " " << endl;

//    cout << "sum1:" << Sum1[0] << endl;

    for (ii=0;ii<chunk_size;ii+=4)
    {

      simd_diff_ += 5*4;

      shift_pointers_2D; //note that these pointers one step further
      load_my_vectors_2D;
 
      calc_zp_2D(C0);
      Sum0 = vec_sldw(B0,C0,1);
      B0=C0;
 
      Sum = vec_add(vec_add(Sum2,Sum1),Sum0);
      vec_st(Sum,0,out);  //commit
      out+=4;
 
      my_0=vec_ld(0,out);
      calc_zz_2D(Sum1);
      calc_zm_2D(C2);
      Sum2 = vec_sldw(B2,C2,3);
      B2 = C2;

//      printf("Sum1:%f %f %f %f\n",Sum1[0],Sum1[1],Sum1[2],Sum1[3]);
//      printf("double_vec:%f %f %f %f\n",double_vec[0],double_vec[1],double_vec[2],double_vec[3]);
//      printf("double_vec2:%f %f %f %f\n",*((double*)(simd_diff_+4*5)),*((double*)(simd_diff_+4*5+2)),*((double*)(simd_diff_+4*5+4)),*((double*)(simd_diff_+4*5+6)));
    }
    out0 += out_jump;
  }
}

void
FGRDiffusionOverlap::compareVoltage(Array3d<double>& A, Array3d<double>& B)
{
   for(int x=0;x<VmBlock_.nx();x++)
   for(int y=0;y<VmBlock_.ny();y++)
   for(int z=0;z<VmBlock_.nz();z++)
   {
     if(A(x,y,z) != B(x,y,z)) printf("comp:%d,%d,%d differ, %g,%g\n",x,y,z,A(x,y,z),B(x,y,z)); 
   }
}


void
FGRDiffusionOverlap::copySlabToBlock()
{
   for(int x=0;x<VmBlock_.nx();x++)
   for(int y=0;y<VmBlock_.ny();y++)
   for(int z=0;z<VmBlock_.nz();z++)
   {
      if ( x == 0)                    VmBlock_(x,y,z)=VmSlab_[0](0,y,z);
      else if( x == VmBlock_.nx()-1 ) VmBlock_(x,y,z)=VmSlab_[1](0,y,z);
      else if( y == 0 )               VmBlock_(x,y,z)=VmSlab_[2](0,x,z);
      else if( y == VmBlock_.ny()-1 ) VmBlock_(x,y,z)=VmSlab_[3](0,x,z);
      else if( z == 0 )               VmBlock_(x,y,z)=VmSlab_[4](0,x,y);
      else if( z == VmBlock_.nz()-1 ) VmBlock_(x,y,z)=VmSlab_[5](0,x,y);
   }
}


void
FGRDiffusionOverlap::clearBlockBD()
{
   for(int x=0;x<VmBlock_.nx();x++)
   for(int y=0;y<VmBlock_.ny();y++)
   for(int z=0;z<VmBlock_.nz();z++)
   {
      if ( x == 0)                    VmBlock_(x,y,z)=0;
      else if( x == VmBlock_.nx()-1 ) VmBlock_(x,y,z)=0;
      else if( y == 0 )               VmBlock_(x,y,z)=0;
      else if( y == VmBlock_.ny()-1 ) VmBlock_(x,y,z)=0;
      else if( z == 0 )               VmBlock_(x,y,z)=0;
      else if( z == VmBlock_.nz()-1 ) VmBlock_(x,y,z)=0;
   }
}
