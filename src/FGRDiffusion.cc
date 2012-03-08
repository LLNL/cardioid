#include "FGRDiffusion.hh"
#include <cassert>
#include "DiffusionUtils.hh"
#include "Anatomy.hh"
#include "Vector.hh"
#include <algorithm>
#include <cstdio>
#include "simd_op.h"

#define check_same
using namespace std;

enum NodeLocation
{ ZZZ, MZZ, ZMZ, PZZ, ZPZ, ZZM, ZZP, MMZ, PMZ, PPZ,
  MPZ, MZM, ZMM, PZM, ZPM, MZP, ZMP, PZP, ZPP
};

namespace
{
   void setGradientWeights(double* grad, int* tissue,
                           NodeLocation n3, NodeLocation n2, NodeLocation n1,
                           NodeLocation n4, NodeLocation n5, NodeLocation n0);
   void f2(unsigned iFace, int tissue[19], double gradPhi[3][19]);
   int isTissue(int cellType){return cellType > 0;}
}

FGRDiffusion::FGRDiffusion(const FGRDiffusionParms& parms,
                           const Anatomy& anatomy)
: localGrid_(DiffusionUtils::findBoundingBox(anatomy)),
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
   
}




void FGRDiffusion::calc(const vector<double>& Vm, vector<double>& dVm, double *recv_buf, int nLocal)
{
   updateVoltageBlock(Vm, recv_buf, nLocal);
   
   int nCell = dVm.size();
   for (int iCell=0; iCell<nCell; ++iCell)
   {
      dVm[iCell] = 0.0;
      int ib = blockIndex_[iCell];
      
      double* phi = & (VmBlock_(ib));
      const double* A = weight_(ib).A;
      for (unsigned ii=0; ii<19; ++ii)
         dVm[iCell] += A[ii] * ( *(phi+offset_[ii]));
      
      dVm[iCell] *= diffusionScale_;
   }
}

void FGRDiffusion::calc_simd(const vector<double>& Vm, vector<double>& dVm, double *recv_buf_, int nLocal)
{
   updateVoltageBlock(Vm, recv_buf_, nLocal);
   int nCell = dVm.size();

   #ifdef check_same
   calc(Vm,dVm,recv_buf_,nLocal); //updateVoltageBlock twice may be ok
   for(int ii=0;ii<20;ii++) std::cout << dVm[ii] << " " ;
   std::cout << std::endl;
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
//       if ((*VmTmp)(ii,jj,kk) !=0 ) cout << VmBlock_(ii,jj,kk) << " ";
     }
     cout << endl;
   }

   Array3d<double> tmp_dVm(VmTmp->nx(),VmTmp->ny(),VmTmp->nz(),0.0);

   uint32_t start = VmTmp->tupleToIndex(1,1,0);
   uint32_t end = VmTmp->tupleToIndex(VmTmp->nx()-2,VmTmp->ny()-2,VmTmp->nz());

//   cout << "start=" << start << endl;
//   cout << "end=" << end << endl;
//   cout << "tmp_dVm.cBlock()=" << tmp_dVm.cBlock() << endl;
//   cout << "tmp_dVm.size=" << tmp_dVm.size() << " " << tmp_dVm.nx() << " " <<  tmp_dVm.ny()  << " " <<  tmp_dVm.nz() <<endl;
//   cout << "VmTmp.size=" << VmTmp->size() << endl;
//   cout << "diffCoef_T2.size=" << diffCoefT2_.size() << " " << diffCoefT2_.nx() << " " << diffCoefT2_.ny() <<  " " << diffCoefT2_.nz() << endl;

//   VmTmp->dump(2,2,10);

 //  diffCoefT2_.dump(2,2,40);

   FGRDiff_simd(start,end-start,VmTmp,tmp_dVm.cBlock());

   if(VmBlock_.nz()%4 != 0) delete VmTmp;

//   tmp_dVm.dump(2,2,10);

   #ifdef check_same
   cout << "checking discrepancy... ";
   for (int ii=0; ii<nCell; ++ii)
   {
      double tmp = dVm[ii];
      dVm[ii] = tmp_dVm(localTuple_[ii].x(),localTuple_[ii].y(),localTuple_[ii].z());
      dVm[ii] *= diffusionScale_;
      if( fabs(tmp - dVm[ii]) > 0.0000001 ) cout << ii << ":" << tmp-dVm[ii] << " " ;
   }
   cout << "Done" << endl;
   for(int ii=0;ii<20;ii++) std::cout << dVm[ii] << " " ;
   std::cout << std::endl;
   #else
   for (int ii=0; ii<nCell; ++ii)
   {
      dVm[ii] = tmp_dVm(localTuple_[ii].x(),localTuple_[ii].y(),localTuple_[ii].z());
      dVm[ii] *= diffusionScale_;
   }
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
   }
   assert(nLocal <= Vm.size());
   for (unsigned ii=nLocal; ii<Vm.size(); ++ii)
   {
      int index = blockIndex_[ii];
      VmBlock_(index) = recv_buf[ii-nLocal];
   }
}

namespace
{
   void setGradientWeights(double* grad, int* tissue,
                           NodeLocation n3, NodeLocation n2, NodeLocation n1,
                           NodeLocation n4, NodeLocation n5, NodeLocation n0)
   {
      int caseNum = tissue[n4] + 2*tissue[n3] + 4*tissue[n0] + 8*tissue[n1];

      double& w0 = grad[n0];
      double& w1 = grad[n1];
      double& w2 = grad[n2];
      double& w3 = grad[n3];
      double& w4 = grad[n4];
      double& w5 = grad[n5];

      switch (caseNum)
      {
        case 0:
         // all weights are zero
         break;
        case 1:
         w5 = 1.0; w4 = -1.0;
         break;
        case 2:
         w2 = 1.0; w3 = -1.0;
         break;
        case 3:
         w2 = w5 = 0.5; w3 = w4 = -0.5;
         break;
        case 4:
         w0 = 1.0; w5 = -1.0;
         break;
        case 5:
         w0 = 0.5; w4 = -0.5;
         break;
        case 6:
         w0 = w2 = 0.5; w3 = w5 = -0.5;
         break;
        case 7:
         w0 = 0.5; w2 = 0.25; w3 = w4 = w5 = -0.25;
         break;
        case 8:
         w1 = 1.0; w2 = -1.0;
         break;
        case 9:
         w1 = w5 = 0.5; w2 = w4 = -0.5;
         break;
        case 10:
         w1 = 0.5; w3 = -0.5;
         break;
        case 11:
         w1 = 0.5; w5 = 0.25; w2 = w3 = w4 = -0.25;
         break;
        case 12:
         w0 = w1 = 0.5; w2 = w5 = -0.5;
         break;
        case 13:
         w4 = -0.5; w2 = -0.25; w0 = w1 = w5 = 0.25;
         break;
        case 14:
         w3 = -0.5; w5 = -0.25; w0 = w1 = w2 = 0.25;
         break;
        case 15:
         w0 = w1 = 0.25; w3 = w4 = -0.25;
         break;
        default:
         assert(false);
      }
   }
}

namespace
{
   void f2(unsigned iFace, int tissue[19], double gradPhi[3][19])
   {
      switch (iFace)
      {
        case 0:
         gradPhi[0][ZZZ] = 1.0;
         gradPhi[0][MZZ] = -1.0;
         setGradientWeights(gradPhi[1], tissue,
                            MMZ, MZZ, MPZ,
                            ZMZ, ZZZ, ZPZ);
         setGradientWeights(gradPhi[2], tissue,
                            ZZM, ZZZ, ZZP,
                            MZM, MZZ, MZP);
         break;
        case 1:
         setGradientWeights(gradPhi[0], tissue,
                            MZZ, ZZZ, PZZ,
                            MMZ, ZMZ, PMZ);
         gradPhi[1][ZZZ] = 1.0;
         gradPhi[1][ZMZ] = -1.0;
         setGradientWeights(gradPhi[2], tissue,
                            ZMM, ZMZ, ZMP,
                            ZZM, ZZZ, ZZP);
         break;
        case 2:
         gradPhi[0][PZZ] = 1.0;
         gradPhi[0][ZZZ] = -1.0;
         setGradientWeights(gradPhi[1], tissue,
                            ZMZ, ZZZ, ZPZ,
                            PMZ, PZZ, PPZ);
         setGradientWeights(gradPhi[2], tissue,
                            PZM, PZZ, PZP,
                            ZZM, ZZZ, ZZP);
         break;
        case 3:
         setGradientWeights(gradPhi[0], tissue,
                            MPZ, ZPZ, PPZ,
                            MZZ, ZZZ, PZZ);
         gradPhi[1][ZPZ] = 1.0;
         gradPhi[1][ZZZ] = -1.0;
         setGradientWeights(gradPhi[2], tissue,
                            ZZM, ZZZ, ZZP,
                            ZPM, ZPZ, ZPP);
         break;
        case 4:
         setGradientWeights(gradPhi[0], tissue,
                            MZM, ZZM, PZM,
                            MZZ, ZZZ, PZZ);
         setGradientWeights(gradPhi[1], tissue,
                            ZMZ, ZZZ, ZPZ,
                            ZMM, ZZM, ZPM);
         gradPhi[2][ZZZ] = 1.0;
         gradPhi[2][ZZM] = -1.0;
         break;
        case 5:
         setGradientWeights(gradPhi[0], tissue,
                            MZZ, ZZZ, PZZ,
                            MZP, ZZP, PZP);
         setGradientWeights(gradPhi[1], tissue,
                            ZMP, ZZP, ZPP,
                            ZMZ, ZZZ, ZPZ);
         gradPhi[2][ZZP] = 1.0;
         gradPhi[2][ZZZ] = -1.0;
         break;
        default:
         assert(false);
      }
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
    int z4 = (int)(zz/4);

    assert(zz > 0);
    assert(zz < Nz2);

    diffCoefT2_(xx,yy,4*19*z4 + 4* 0 + zz%4 ) =weight_(xx,yy,zz-1).A[ZZP];
    diffCoefT2_(xx,yy,4*19*z4 + 4* 1 + zz%4 ) =weight_(xx,yy,zz-1).A[ZMP];
    diffCoefT2_(xx,yy,4*19*z4 + 4* 2 + zz%4 ) =weight_(xx,yy,zz-1).A[MZP];
    diffCoefT2_(xx,yy,4*19*z4 + 4* 3 + zz%4 ) =weight_(xx,yy,zz-1).A[ZPP];
    diffCoefT2_(xx,yy,4*19*z4 + 4* 4 + zz%4 ) =weight_(xx,yy,zz-1).A[PZP];
                                 
    diffCoefT2_(xx,yy,4*19*z4 + 4* 5 + zz%4 ) =weight_(xx,yy,zz).A[ZZZ];
    diffCoefT2_(xx,yy,4*19*z4 + 4* 6 + zz%4 ) =weight_(xx,yy,zz).A[PMZ];
    diffCoefT2_(xx,yy,4*19*z4 + 4* 7 + zz%4 ) =weight_(xx,yy,zz).A[MMZ];
    diffCoefT2_(xx,yy,4*19*z4 + 4* 8 + zz%4 ) =weight_(xx,yy,zz).A[MPZ];
    diffCoefT2_(xx,yy,4*19*z4 + 4* 9 + zz%4 ) =weight_(xx,yy,zz).A[PPZ];
    diffCoefT2_(xx,yy,4*19*z4 + 4*10 + zz%4 ) =weight_(xx,yy,zz).A[ZMZ];
    diffCoefT2_(xx,yy,4*19*z4 + 4*11 + zz%4 ) =weight_(xx,yy,zz).A[MZZ];
    diffCoefT2_(xx,yy,4*19*z4 + 4*12 + zz%4 ) =weight_(xx,yy,zz).A[ZPZ];
    diffCoefT2_(xx,yy,4*19*z4 + 4*13 + zz%4 ) =weight_(xx,yy,zz).A[PZZ];

    diffCoefT2_(xx,yy,4*19*z4 + 4*14 + zz%4 ) =weight_(xx,yy,zz+1).A[ZZM];
    diffCoefT2_(xx,yy,4*19*z4 + 4*15 + zz%4 ) =weight_(xx,yy,zz+1).A[ZMM];
    diffCoefT2_(xx,yy,4*19*z4 + 4*16 + zz%4 ) =weight_(xx,yy,zz+1).A[MZM];
    diffCoefT2_(xx,yy,4*19*z4 + 4*17 + zz%4 ) =weight_(xx,yy,zz+1).A[ZPM];
    diffCoefT2_(xx,yy,4*19*z4 + 4*18 + zz%4 ) =weight_(xx,yy,zz+1).A[PZM];

//    diffCoefT2_(xx,yy,4*19*z4 + 4* 0 + zz%4 ) =weight_(xx,yy,zz-1).A[ZZP]=rand();
//    diffCoefT2_(xx,yy,4*19*z4 + 4* 1 + zz%4 ) =weight_(xx,yy,zz-1).A[ZMP]=rand();
//    diffCoefT2_(xx,yy,4*19*z4 + 4* 2 + zz%4 ) =weight_(xx,yy,zz-1).A[MZP]=rand();
//    diffCoefT2_(xx,yy,4*19*z4 + 4* 3 + zz%4 ) =weight_(xx,yy,zz-1).A[ZPP]=rand();
//    diffCoefT2_(xx,yy,4*19*z4 + 4* 4 + zz%4 ) =weight_(xx,yy,zz-1).A[PZP]=rand();
//                                 
//    diffCoefT2_(xx,yy,4*19*z4 + 4* 5 + zz%4 ) =weight_(xx,yy,zz).A[ZZZ]=rand();
//    diffCoefT2_(xx,yy,4*19*z4 + 4* 6 + zz%4 ) =weight_(xx,yy,zz).A[PMZ]=rand();
//    diffCoefT2_(xx,yy,4*19*z4 + 4* 7 + zz%4 ) =weight_(xx,yy,zz).A[MMZ]=rand();
//    diffCoefT2_(xx,yy,4*19*z4 + 4* 8 + zz%4 ) =weight_(xx,yy,zz).A[MPZ]=rand();
//    diffCoefT2_(xx,yy,4*19*z4 + 4* 9 + zz%4 ) =weight_(xx,yy,zz).A[PPZ]=rand();
//    diffCoefT2_(xx,yy,4*19*z4 + 4*10 + zz%4 ) =weight_(xx,yy,zz).A[ZMZ]=rand();
//    diffCoefT2_(xx,yy,4*19*z4 + 4*11 + zz%4 ) =weight_(xx,yy,zz).A[MZZ]=rand();
//    diffCoefT2_(xx,yy,4*19*z4 + 4*12 + zz%4 ) =weight_(xx,yy,zz).A[ZPZ]=rand();
//    diffCoefT2_(xx,yy,4*19*z4 + 4*13 + zz%4 ) =weight_(xx,yy,zz).A[PZZ]=rand();
//
//    diffCoefT2_(xx,yy,4*19*z4 + 4*14 + zz%4 ) =weight_(xx,yy,zz+1).A[ZZM]=rand();
//    diffCoefT2_(xx,yy,4*19*z4 + 4*15 + zz%4 ) =weight_(xx,yy,zz+1).A[ZMM]=rand();
//    diffCoefT2_(xx,yy,4*19*z4 + 4*16 + zz%4 ) =weight_(xx,yy,zz+1).A[MZM]=rand();
//    diffCoefT2_(xx,yy,4*19*z4 + 4*17 + zz%4 ) =weight_(xx,yy,zz+1).A[ZPM]=rand();
//    diffCoefT2_(xx,yy,4*19*z4 + 4*18 + zz%4 ) =weight_(xx,yy,zz+1).A[PZM]=rand();

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

  const double* VmM = VmTmp->cBlock();

  int xm1ym1z_ = ((-1) *Ny2 + (-1)) * Nz2;
  int xm1yz_ =   ((-1) *Ny2 + ( 0)) * Nz2;
  int xm1yp1z_ = ((-1) *Ny2 + (+1)) * Nz2;
  int xym1z_ =   ((0 ) *Ny2 + (-1)) * Nz2;
  int xyp1z_ =   ((0 ) *Ny2 + (+1)) * Nz2;
  int xp1ym1z_ = ((+1) *Ny2 + (-1)) * Nz2;
  int xp1yz_ =   ((+1) *Ny2 + ( 0)) * Nz2;
  int xp1yp1z_ = ((+1) *Ny2 + (+1)) * Nz2;

  const double* phi_xm1_ym1_z = VmM + start + xm1ym1z_;
  const double* phi_xm1_y_z   = VmM + start + xm1yz_;
  const double* phi_xm1_yp1_z = VmM + start + xm1yp1z_;
  const double* phi_x_ym1_z   = VmM + start + xym1z_;
  const double* phi_x_y_z     = VmM + start ;
  const double* phi_x_yp1_z   = VmM + start + xyp1z_;
  const double* phi_xp1_ym1_z = VmM + start + xp1ym1z_;
  const double* phi_xp1_y_z   = VmM + start + xp1yz_;
  const double* phi_xp1_yp1_z = VmM + start + xp1yp1z_;

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
 
  double *simd_diff_ = diffCoefT2_.cBlock() + start*19 - 4;
//  printf("simd_diff:%10f %10f %10f %10f\n",*(simd_diff_+4),*(simd_diff_+5),*(simd_diff_+6),*(simd_diff_+7));
//  printf("simd_diff:%10f %10f %10f %10f\n",*(simd_diff_+8),*(simd_diff_+9),*(simd_diff_+10),*(simd_diff_+11));

  B0 =  vec_madd(vec_ld(0,simd_diff_+=4)  , my_xp1_y_z,
        vec_madd(vec_ld(0,simd_diff_+=4)  , my_x_yp1_z,
        vec_madd(vec_ld(0,simd_diff_+=4)  , my_xm1_y_z,
        vec_madd(vec_ld(0,simd_diff_+=4)  , my_x_ym1_z,
        vec_mul( vec_ld(0,simd_diff_+=4)  , my_x_y_z)))));

  B1 = vec_madd( vec_ld(0,simd_diff_+=4)  , my_xp1_y_z,
       vec_madd( vec_ld(0,simd_diff_+=4)  , my_x_yp1_z,
       vec_madd( vec_ld(0,simd_diff_+=4)  , my_xm1_y_z,
       vec_madd( vec_ld(0,simd_diff_+=4)  , my_x_ym1_z,
       vec_madd( vec_ld(0,simd_diff_+=4)  , my_xp1_yp1_z,
       vec_madd( vec_ld(0,simd_diff_+=4)  , my_xm1_yp1_z,
       vec_madd( vec_ld(0,simd_diff_+=4)  , my_xm1_ym1_z,
       vec_madd( vec_ld(0,simd_diff_+=4)  , my_xp1_ym1_z,
       vec_mul ( vec_ld(0,simd_diff_+=4)  , my_x_y_z)))))))));
//  printf("%f %f %f %f\n",my_x_y_z.v[0], my_x_y_z.v[1], my_x_y_z.v[2], my_x_y_z.v[3]);
//  printf("%f %f %f %f\n",B1.v[0],B1.v[1],B1.v[2],B1.v[3]);

  B2 = vec_madd( vec_ld(0,simd_diff_+=4)  , my_xp1_y_z,
       vec_madd( vec_ld(0,simd_diff_+=4)  , my_x_yp1_z,
       vec_madd( vec_ld(0,simd_diff_+=4)  , my_xm1_y_z,
       vec_madd( vec_ld(0,simd_diff_+=4)  , my_x_ym1_z,
       vec_mul ( vec_ld(0,simd_diff_+=4)  , my_x_y_z)))));

  Sum2 = vec_sldw(my_zero_vec,B2,3);

  for(ii=1;ii<chunk_size;ii+=4)  //start must be that os x chunk.
  {
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


  C0 =  vec_madd(vec_ld(0,simd_diff_+=4)  , my_xp1_y_z,
        vec_madd(vec_ld(0,simd_diff_+=4)  , my_x_yp1_z,
        vec_madd(vec_ld(0,simd_diff_+=4)  , my_xm1_y_z,
        vec_madd(vec_ld(0,simd_diff_+=4)  , my_x_ym1_z,
        vec_mul( vec_ld(0,simd_diff_+=4)  , my_x_y_z)))));

    Sum0 = vec_sldw(B0,C0,1);
    B0 = C0;

//  printf("%f %f %f %f\n",Sum0.v[0],Sum0.v[1],Sum0.v[2],Sum0.v[3]);
//  printf("%f %f %f %f\n",Sum2.v[0],Sum2.v[1],Sum2.v[2],Sum2.v[3]);
//  printf("%f %f %f %f\n",(*simd_diff_),*(simd_diff_+1),*(simd_diff_+2),*(simd_diff_+3));
    // commit the previous
    Sum = vec_add(vec_add(Sum0,B1),Sum2);

    vec_st(Sum,0,out);
//  printf("out:%f %f %f %f\n",out[0],out[1],out[2],out[3]);

  B1 = vec_madd( vec_ld(0,simd_diff_+=4)  , my_xp1_y_z,
       vec_madd( vec_ld(0,simd_diff_+=4)  , my_x_yp1_z,
       vec_madd( vec_ld(0,simd_diff_+=4)  , my_xm1_y_z,
       vec_madd( vec_ld(0,simd_diff_+=4)  , my_x_ym1_z,
       vec_madd( vec_ld(0,simd_diff_+=4)  , my_xp1_yp1_z,
       vec_madd( vec_ld(0,simd_diff_+=4)  , my_xm1_yp1_z,
       vec_madd( vec_ld(0,simd_diff_+=4)  , my_xm1_ym1_z,
       vec_madd( vec_ld(0,simd_diff_+=4)  , my_xp1_ym1_z,
       vec_mul ( vec_ld(0,simd_diff_+=4)  , my_x_y_z)))))))));

//  printf("B1 :%f %f %f %f\n",B1.v[0],B1.v[1],B1.v[2],B1.v[3]);
//  printf("C0 :%f %f %f %f\n",C0.v[0],C0.v[1],C0.v[2],C0.v[3]);
//  printf("C2 :%f %f %f %f\n",C2.v[0],C2.v[1],C2.v[2],C2.v[3]);
//  printf("xyz:%f %f %f %f\n",my_x_y_z.v[0], my_x_y_z.v[1], my_x_y_z.v[2], my_x_y_z.v[3]);


  C2 = vec_madd( vec_ld(0,simd_diff_+=4)  , my_xp1_y_z,
       vec_madd( vec_ld(0,simd_diff_+=4)  , my_x_yp1_z,
       vec_madd( vec_ld(0,simd_diff_+=4)  , my_xm1_y_z,
       vec_madd( vec_ld(0,simd_diff_+=4)  , my_x_ym1_z,
       vec_mul ( vec_ld(0,simd_diff_+=4)  , my_x_y_z)))));

    Sum2 = vec_sldw(B2,C2,3);
    B2 = C2;

   out +=4;

  }
}

