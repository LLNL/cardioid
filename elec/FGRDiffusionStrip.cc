#include "FGRDiffusionStrip.hh"
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


FGRDiffusionStrip::FGRDiffusionStrip(const FGRDiffusionParms& parms,
                           const Anatomy& anatomy,
                           const ThreadTeam& threadInfo,
                           const ThreadTeam& reactionThreadInfo)
: Diffusion(parms.diffusionScale_),
  nLocal_(anatomy.nLocal()),
  localGrid_(DiffusionUtils::findBoundingBox_simd(anatomy, parms.printBBox_)),
  threadInfo_(threadInfo),
  reactionThreadInfo_(reactionThreadInfo)
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

   mkOffsets(threadOffset_,     anatomy.nLocal(),  threadInfo_);
   mkOffsets(localCopyOffset_,  anatomy.nLocal(),  reactionThreadInfo_);
   mkOffsets(remoteCopyOffset_, anatomy.nRemote(), threadInfo_);

   //simd thread offsets
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

   assert(nz%4 == 0);
   dVmBlock_.resize(nx,ny,nz + (nz%4==0 ? 0:4-(nz%4)));

#define STRIP_DIFFUSION_DEBUG 0

   if(STRIP_DIFFUSION_DEBUG)
     fprintf(stderr,"%s:%06d: Making compute strips\n",__FILE__,__LINE__);
   /* Make compute strips */ {
     int nt = threadInfo.nThreads();
     int *board = new int[nx*ny];
     int ctot = 0,c0,c1;
     int is0,is1;

     if(STRIP_DIFFUSION_DEBUG)
       fprintf(stderr,"%s:%d: nt = %d, nx = %d, ny = %d, nx*ny = %d\n",
	       __FILE__,__LINE__,nt,nx,ny,nx*ny);

     strip_begin.resize(nt);
     strip_end.resize(nt);

     /* Initialize board (xy-table of emptiness of z-columns) */ {
       const vector<AnatomyCell>& cell = anatomy.cellArray();

       for(int i = 0; i<nx*ny; i++) board[i] = 0;

       for (unsigned ii=0; ii<anatomy.nLocal(); ++ii) {
	 Tuple globalTuple = anatomy.globalTuple(ii);
	 Tuple ll = localGrid_.localTuple(globalTuple);

	 if(ll.x() < 0 || ll.x() >= nx ||
	    ll.y() < 0 || ll.y() >= ny) {
	   fprintf(stderr,"%s:%d: error, x or y out of range: x=%d, y=%d, nx=%d, ny=%d\n",
		   __FILE__,__LINE__,(int) ll.x(),(int) ll.y(),nx,ny);
	   exit(1);
	 }

	 if(ll.x() < 1 || ll.x() >= nx-1 ||
	    ll.y() < 1 || ll.y() >= ny-1) {
	   fprintf(stderr,"%s:%d: Actually, I expected all local cells to be in the interior of [0,nx-1]x[0,ny-1]  x=%d, y=%d, nx=%d, ny=%d\n",
		   __FILE__,__LINE__,(int) ll.x(),(int) ll.y(),nx,ny);
	   exit(1);
	 }

	 board[ll.y() + ny*ll.x()] = 1;
       }
     }
     for(int i = 0; i<nx*ny; i++)
       if(board[i] == 1) ctot = ctot + 1;

     if(STRIP_DIFFUSION_DEBUG) {
       fprintf(stderr,"    : ");
       for(int j = 0; j<ny; j++)
	 fprintf(stderr,"%2d",j);
       fprintf(stderr,"\n");
       for(int i = 0; i<nx; i++) {
	 fprintf(stderr,"%4d: ",i);
	 for(int j = 0; j<ny; j++)
	   fprintf(stderr,"%2d",board[j+i*ny]);
	 fprintf(stderr,"\n");
       }
     }

     if(STRIP_DIFFUSION_DEBUG)
       fprintf(stderr,"%s:%d: Board initialized, ctot = %d\n",__FILE__,__LINE__,ctot);

     for(int tid = 0; tid<nt; tid++) {
       
       /* Split non-empty z-columns over threads */ {
	 int q = ctot / nt,r = ctot % nt;
	 if(tid < r) {
	   is0 = (q+1)*tid;
	   is1 = is0 + q + 1;
	 } else {
	   is0 = q*tid + r;
	   is1 = is0 + q;
	 }
       }
       
       if(is1 > is0) {
	 /* Find start and end column for this thread */ {
	   int c = 0;
	   for(int i = 0; i<nx*ny; i++) {
	     if(board[i] > 0) {
	       if(c == is0) c0 = i;
	       c = c + 1;
	       if(c == is1) c1 = i;
	     }
	   }
	 }

	 if(STRIP_DIFFUSION_DEBUG)
	   fprintf(stderr,"%s:%d: tid = %2d, is0=%d, is1=%d, c0=%d, c1=%d\n",__FILE__,__LINE__,
		   tid,is0,is1,c0,c1);
	 /* Build list of chunks (conscutive non-empty columns) */ {
	   int empty = 0;
	   strip_begin[tid].push_back(c0);
	   for(int i = c0; i<=c1; i++) {
	     if(empty == 1) {
	       if(board[i] == 1) {
		 strip_begin[tid].push_back(i);
		 empty = 0;
	       }
	     } else /* if(empty == 0) */ {
	       if(board[i] == 0) {
		 strip_end[tid].push_back(i-1);
		 empty = 1;
		 if(STRIP_DIFFUSION_DEBUG)
		   fprintf(stderr,"    Strip %d: begin=%d  end=%d\n",
			   (int) strip_end[tid].size(),
			   
			   strip_begin[tid][strip_end[tid].size()-1],
			   strip_end[tid][strip_end[tid].size()-1]);
	       }
	     }
	   }
	   if(strip_end[tid].size() < strip_begin[tid].size()) {
	     strip_end[tid].push_back(c1);
	     if(STRIP_DIFFUSION_DEBUG)
	       fprintf(stderr,"    Strip %d: begin=%d  end=%d\n",
		       (int) strip_end[tid].size(),strip_begin[tid][strip_end[tid].size()-1],c1);
	   }
	 }
       }
       if(STRIP_DIFFUSION_DEBUG)
	 fprintf(stderr,"Strip construction completed.\n\n");
     }
     /* Test strip distribution for overlap and completeness */ {
       if(STRIP_DIFFUSION_DEBUG) fprintf(stderr,"nx = %d,  ny = %d, nx*ny = %d\n",nx,ny,nx*ny);
       int *tag = new int[nx*ny];

       for(int i = 0; i<nx*ny; i++) tag[i] = 0;
       if(STRIP_DIFFUSION_DEBUG)
	 fprintf(stderr,"tag allocated and initialized. nx=%d ny=%d nx*ny=%d\n",nx,ny,nx*ny);

       for(int tid = 0; tid<nt; tid++) {
	 for(int is = 0; is<(int) strip_begin[tid].size(); is++) {
	   int j0 = strip_begin[tid][is];
	   int j1 = strip_end[tid][is];
	   if(j0 >= 0 && j1 >= j0 && j1 < nx*ny)
	     for(int j = j0; j<=j1; j++)
	       tag[j] = tag[j] + 1;
	   else {
	     fprintf(stderr,
		     "%s:%d: error, strip index out of range.\n"
		     "nx=%d ny=%d nx*ny=%d tid=%d nt=%d is=%d j0=%d j1=%d\n",
		     __FILE__,__LINE__,
		     nx,ny,nx*ny,tid,nt,is,j0,j1);
	     exit(1);
	   }
	 }
       }

       for(int i = 0; i<nx*ny; i++)
	 if(tag[i] != board[i]) {
	   fprintf(stderr,
		   "%s:%d: error, tag count is wrong. nx=%d ny=%d nx*ny=%d i=%d tag=%d board=%d\n",
		   __FILE__,__LINE__,
		   nx,ny,nx*ny,i,tag[i],board[i]);
	   exit(1);
	 }

       delete[] tag;
     }

     delete[] board;
   }


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

void FGRDiffusionStrip::updateLocalVoltage(ro_mgarray_ptr<double> VmLocal_managed)
{
   startTimer(FGR_ArrayLocal2MatrixTimer);
   ro_array_ptr<double> VmLocal = VmLocal_managed.useOn(CPU);
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

void FGRDiffusionStrip::updateRemoteVoltage(ro_mgarray_ptr<double> VmRemote_managed)
{
   startTimer(FGR_ArrayRemote2MatrixTimer);
   ro_array_ptr<double> VmRemote = VmRemote_managed.useOn(CPU);
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

/** threaded simd version */
void FGRDiffusionStrip::calc(rw_mgarray_ptr<double> dVm_managed)
{
   rw_array_ptr<double> dVm = dVm_managed.useOn(CPU);
   int tid = threadInfo_.teamRank();
   startTimer(FGR_AlignCopyTimer);

   Array3d<double> *VmTmp = &(VmBlock_);

   //make sure z is multiple of 4
   if(VmBlock_.nz()%4 != 0)
   {
      assert(false); // this shouldn't ever happen.
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
   stopTimer(FGR_AlignCopyTimer);
   
   startTimer(FGR_StencilTimer);

   //uint32_t begin = VmTmp->tupleToIndex(threadOffsetSimd_[tid],1,0);
   //uint32_t end = VmTmp->tupleToIndex(threadOffsetSimd_[tid+1]-1,VmTmp->ny()-2,VmTmp->nz());
   //   printf("simd version:%d-%d\n",begin,end);
   //if (threadOffsetSimd_[tid] < threadOffsetSimd_[tid+1] )
   //   FGRDiff_simd_thread(begin,end-begin,VmTmp,dVmBlock_.cBlock());

   //if (threadOffsetSimd_[tid] < threadOffsetSimd_[tid+1] )
   FGRDiff_simd_thread(threadOffsetSimd_[tid] ,
		       threadOffsetSimd_[tid+1],
		       VmTmp,dVmBlock_.cBlock());

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
void FGRDiffusionStrip::buildTupleArray(const Anatomy& anatomy)
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

void FGRDiffusionStrip::buildBlockIndex(const Anatomy& anatomy)
{
   blockIndex_.resize(anatomy.size());
   for (unsigned ii=0; ii<anatomy.size(); ++ii)
   {
      Tuple globalTuple = anatomy.globalTuple(ii);
      Tuple ll = localGrid_.localTuple(globalTuple);
      blockIndex_[ii] = VmBlock_.tupleToIndex(ll.x(), ll.y(), ll.z());
   }
}

//void FGRDiffusionStrip::precomputeCoefficientsPositive(const Anatomy& anatomy)
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


void FGRDiffusionStrip::precomputeCoefficients(const Anatomy& anatomy)
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
//   printAllWeights(tissueBlk);
}

void FGRDiffusionStrip::mkTissueArray(
   const Array3d<int>& tissueBlk, int ib, int* tissue)
{
   for (unsigned ii=0; ii<19; ++ii)
      tissue[ii] = tissueBlk(ib + offset_[ii]);
}


Vector FGRDiffusionStrip::f1(int ib, int iFace, const Vector& h,
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


void FGRDiffusionStrip::printAllWeights(const Array3d<int>& tissue)
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

void FGRDiffusionStrip::reorder_Coeff()
{
  uint32_t idx=0,ll,ii,jj,kk;
  const uint32_t Nx2 = localGrid_.nx();
  const uint32_t Ny2 = localGrid_.ny();
  const uint32_t Nz2 = localGrid_.nz() + ((localGrid_.nz()%4)==0 ? 0 : (4-localGrid_.nz()%4));

  diffCoefT2_.resize(Nx2,Ny2,Nz2*20,0.0);
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

    diffCoefT2_(xx,yy,4*20*zp4 + 4* 0 + zp%4 ) =weight_(xx,yy,zz).A[ZZP];
    diffCoefT2_(xx,yy,4*20*zp4 + 4* 1 + zp%4 ) =weight_(xx,yy,zz).A[ZMP];
    diffCoefT2_(xx,yy,4*20*zp4 + 4* 2 + zp%4 ) =weight_(xx,yy,zz).A[MZP];
    diffCoefT2_(xx,yy,4*20*zp4 + 4* 3 + zp%4 ) =weight_(xx,yy,zz).A[ZPP];
    diffCoefT2_(xx,yy,4*20*zp4 + 4* 6 + zp%4 ) =weight_(xx,yy,zz).A[PZP];
                                 
    *((double*)&(diffCoefT2_(xx,yy,4*20*zz4 + 4* 4 + (zz%4)*2))) = A0_(xx,yy,zz);
//    diffCoefT2_(xx,yy,4*20*zz4 + 4* 5 + zz%4 ) =weight_(xx,yy,zz).A[ZZZ];
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

  }
}



// simdiazed version
// 'start' cannot be in the middle. 'start' must point to (n,m,0). 
void
FGRDiffusionStrip::FGRDiff_simd_thread(const uint32_t bx,const int32_t ex, Array3d<double>* VmTmp, double* out0)
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

  /* for(int xx=bx;xx<ex;xx++) */
  int tid = threadInfo_.teamRank();
  std::vector<int> &s_begin = strip_begin[tid],&s_end = strip_end[tid];
  for(int istrip = 0; istrip < s_begin.size(); istrip++)
  {
    /*
      int start = VmTmp->tupleToIndex(xx,1,0);
      int end   = VmTmp->tupleToIndex(xx,Ny2-2,Nz2);
      int chunk_size = end - start;
    */
    int idx0 = s_begin[istrip], idx1 = s_end[istrip];
    int ix0 = idx0 / Ny2, iy0 = idx0 % Ny2;
    int ix1 = idx1 / Ny2, iy1 = idx1 % Ny2;
    int start = VmTmp->tupleToIndex(ix0,iy0,0);
    int end   = VmTmp->tupleToIndex(ix1,iy1,0);
    int chunk_size = end - start + Nz2;

    if(ix0 != ix1) {
      fprintf(stderr,
	      "%s:%d: Gddmmt!"
	      "There is supposed to be empty columns at y=0 and y=ny-1.\n",
	      __FILE__,__LINE__);
      assert(ix0 == ix1);
    }

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
}


