#include "FGRDiffusion.hh"
#include <cassert>
#include "DiffusionUtils.hh"
#include "Anatomy.hh"
#include "Vector.hh"
#include <algorithm>

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
   
}




void FGRDiffusion::calc(const vector<double>& Vm, vector<double>& dVm)
{
   updateVoltageBlock(Vm);
   
   for (unsigned iCell=0; iCell<Vm.size(); ++iCell)
   {
      dVm[iCell] = 0.0;
      int ib = blockIndex_[iCell];
      
      double* phi = & (VmBlock_(ib));
      const double* A = weight_(ib).A;
      for (unsigned ii=0; ii<19; ++ii)
         dVm[iCell] -= A[ii] * ( *(phi+offset_[ii]));
      
      dVm[iCell] *= diffusionScale_;
   }
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

         Vector sigmaTimesS = f1(ib, iFace, hInv, sigmaBlk);
         double gradPhi[3][19] = {0};
         f2(iFace, tissue, gradPhi);
         
         for (unsigned ii=0; ii<19; ++ii)
            for (unsigned jj=0; jj<3; ++jj)
               weight_(ib).A[ii] += sigmaTimesS[jj] * gradPhi[jj][ii] * hInv[jj];
      }
   }
//   printAllWeights(tissueBlk);
}


void FGRDiffusion::updateVoltageBlock(const vector<double>& Vm)
{
   for (unsigned ii=0; ii<Vm.size(); ++ii)
   {
      int index = blockIndex_[ii];
      VmBlock_(index) = Vm[ii];
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
         w1 = 1.0; w3 = -1.0;
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
         w5 = -0.5; w5 = -0.25; w0 = w1 = w2 = 0.25;
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
                            ZMZ, ZZZ, ZPZ,
                            PMZ, PZZ, PPZ);
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
                            MMZ, MZZ, MPZ,
                            ZMZ, ZZZ, ZPZ);
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


Vector FGRDiffusion::f1(int ib, int iFace, const Vector& hInv,
                        const Array3d<SymmetricTensor>& sigmaBlk)
{
   SymmetricTensor
      sigma = (sigmaBlk(ib) + sigmaBlk(ib+faceNbrOffset_[iFace])) / 2.0;
   Vector S(0, 0, 0);
   switch (iFace)
   {
     case 0:
      S[0] = -hInv[1]*hInv[2];
      break;
     case 1:
      S[1] = -hInv[0]*hInv[2];
      break;
     case 2:
      S[0] = hInv[1]*hInv[2];
      break;
     case 3:
      S[1] = hInv[0]*hInv[2];
      break;
     case 4:
      S[2] = -hInv[0]*hInv[1];
      break;
     case 5:
      S[2] = hInv[0]*hInv[1];
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
