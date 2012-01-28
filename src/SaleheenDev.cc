#include "SaleheenDev.hh"

#include "Anatomy.hh"
#include "SymmetricTensor.hh"
#include <algorithm>
#include <iostream>

using namespace std;

/** We want to find the boundingBox such that any stencil point of any
 *  local atom is in the box.  It is not sufficient merely to iterate all
 *  of the local and remote atoms and find the maximum extent.  There
 *  may be local cells that are on the outer or inner walls of the
 *  heart.  Such cells will have no remote cells to satisfy their
 *  stencil.  Therefore, the safe bet is to iterate the local cells and
 *  add the stencil size in each direction.
 */
namespace
{
   LocalGrid findBoundingBox(const Anatomy& anatomy)
   {
      assert(anatomy.nLocal() > 0);
      Tuple globalTuple = anatomy.globalTuple(0);
      int xMin = globalTuple.x();
      int yMin = globalTuple.y();
      int zMin = globalTuple.z();
      int xMax = globalTuple.x();
      int yMax = globalTuple.y();
      int zMax = globalTuple.z();
      
      for (unsigned ii=1; ii<anatomy.nLocal(); ++ii)
      {
         Tuple globalTuple = anatomy.globalTuple(ii);
         xMin = min(xMin, globalTuple.x());
         yMin = min(yMin, globalTuple.y());
         zMin = min(zMin, globalTuple.z());
         xMax = max(xMax, globalTuple.x());
         yMax = max(yMax, globalTuple.y());
         zMax = max(zMax, globalTuple.z());
      }
      
      int stencilSize = 1;
      
      int nx = 2*stencilSize + xMax - xMin + 1;
      int ny = 2*stencilSize + yMax - yMin + 1;
      int nz = 2*stencilSize + zMax - zMin + 1;
      xMin -= stencilSize;
      yMin -= stencilSize;
      zMin -= stencilSize;
      
      return LocalGrid(nx, ny, nz, xMin, yMin, zMin);
   }
}

/** Helper function to compute the finite difference approximation of
 *  the gradient of the conductivity.  Uses either a one-sided or
 *  two-sided approximation according to the cell type of the
 *  neighboring points.  In the case that both nbr points are non-tissue
 *  the gradient is taken to be zero.
 */
namespace
{
   double dSigma(double ss, double ssPlus, double ssMinus, double hhInv, int typePlus, int typeMinus)
   {
      int n=0;
      double tmp=0;
      if (typePlus != 0)
      {
         tmp += ssPlus - ss;
         ++n;
      }
      if (typeMinus != 0)
      {
         tmp += ss - ssMinus;
         ++n;
      }
      switch (n)
      {
        case 0:
         break;
        case 1:
         tmp *= hhInv;
         break;
        case 2:
         tmp *= 0.5*hhInv;
         break;
        default:
         assert(false);
      }

      return tmp;
   }
}



SaleheenDev::SaleheenDev(
   const SaleheenDevParms& parms,
   const Anatomy& anatomy)
: localGrid_(findBoundingBox(anatomy)),
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
   
   diffIntra_.resize(nx, ny, nz);
   VmBlock_.resize(nx, ny, nz);

   buildTupleArray(anatomy);
   buildBlockIndex(anatomy);
   precomputeCoefficients(anatomy);

   int base = VmBlock_.tupleToIndex(1, 1, 1);
   xp1yz_   = VmBlock_.tupleToIndex(2, 1, 1) - base;
   xyp1z_   = VmBlock_.tupleToIndex(1, 2, 1) - base;
   xm1yz_   = VmBlock_.tupleToIndex(0, 1, 1) - base;
   xym1z_   = VmBlock_.tupleToIndex(1, 0, 1) - base;
   xp1yp1z_ = VmBlock_.tupleToIndex(2, 2, 1) - base;
   xm1yp1z_ = VmBlock_.tupleToIndex(0, 2, 1) - base;
   xm1ym1z_ = VmBlock_.tupleToIndex(0, 0, 1) - base;
   xp1ym1z_ = VmBlock_.tupleToIndex(2, 0, 1) - base;
   xyzp1_   = VmBlock_.tupleToIndex(1, 1, 2) - base;
   xyzm1_   = VmBlock_.tupleToIndex(1, 1, 0) - base;
   xyp1zp1_ = VmBlock_.tupleToIndex(1, 2, 2) - base;
   xyp1zm1_ = VmBlock_.tupleToIndex(1, 2, 0) - base;
   xym1zm1_ = VmBlock_.tupleToIndex(1, 0, 0) - base;
   xym1zp1_ = VmBlock_.tupleToIndex(1, 0, 2) - base;
   xp1yzp1_ = VmBlock_.tupleToIndex(2, 1, 2) - base;
   xm1yzp1_ = VmBlock_.tupleToIndex(0, 1, 2) - base;
   xm1yzm1_ = VmBlock_.tupleToIndex(0, 1, 0) - base;
   xp1yzm1_ = VmBlock_.tupleToIndex(2, 1, 0) - base;

}

   
void SaleheenDev::calc(
   const vector<double>& Vm, vector<double>& dVm)
{
   updateVoltageBlock(Vm);

   int n = dVm.size();
#pragma omp parallel for
   for (int ii=0; ii<n; ++ii)
   {
      dVm[ii] = boundaryFDLaplacianSaleheen98SumPhi(blockIndex_[ii]);
      dVm[ii] *= diffusionScale_;
   }
}
/** We're building the localTuple array only for local cells.  We can't
 * do stencil operations on remote particles so we shouldn't need
 * tuples.  We can use block indices instead.
 */
void SaleheenDev::buildTupleArray(const Anatomy& anatomy)
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

void SaleheenDev::buildBlockIndex(const Anatomy& anatomy)
{
   blockIndex_.resize(anatomy.size());
   for (unsigned ii=0; ii<anatomy.size(); ++ii)
   {
      Tuple globalTuple = anatomy.globalTuple(ii);
      Tuple ll = localGrid_.localTuple(globalTuple);
      blockIndex_[ii] = VmBlock_.tupleToIndex(ll.x(), ll.y(), ll.z());
   }
}


/** This routine pre-computes the diffusion weights for the finitie
 *  difference approximation given in Saleheen98.
 *
 *  To enforce the zero flux boundary condition for the monodomain
 *  problem we explicitly set all conductivities to zero and then
 *  compute conductivity only for actual tissue cells in the anatomy
 *  cell array.
 *
 *  We also ensure that the approximation for the spatial gradient of
 *  the conductivities uses a one-sided difference for cells on the
 *  organ surface.
 */
void
SaleheenDev::precomputeCoefficients(const Anatomy& anatomy)
{
   unsigned nx = localGrid_.nx();
   unsigned ny = localGrid_.ny();
   unsigned nz = localGrid_.nz();
   unsigned nxGlobal = anatomy.nx();
   unsigned nyGlobal = anatomy.ny();
   unsigned nzGlobal = anatomy.nz();

   SymmetricTensor sigmaZero = {0};
   Array3d<SymmetricTensor> sigmaMintra(nx, ny, nz, sigmaZero);
   Array3d<int> tissue(nx, ny, nz, 0);

   const vector<AnatomyCell>& cell = anatomy.cellArray();
   for (unsigned ii=0; ii<anatomy.size(); ++ii)
   {
      unsigned ib = blockIndex_[ii];
      sigmaMintra(ib) = anatomy.conductivity(ii);
      tissue(ib) = anatomy.cellType(ii);
   }

   double dxInv = 1.0/anatomy.dx();
   double dyInv = 1.0/anatomy.dy();
   double dzInv = 1.0/anatomy.dz();

   for (unsigned ii=0; ii<anatomy.nLocal(); ++ii)
   {
      unsigned ix = localTuple_[ii].x();
      unsigned iy = localTuple_[ii].y();
      unsigned iz = localTuple_[ii].z();
      
      // compute diffIntra_(ix, iy, iz)
      const int*** tt = (const int***) tissue.cArray();
      const SymmetricTensor*** ss =
         (const SymmetricTensor***) sigmaMintra.cArray();
      boundaryFDLaplacianSaleheen98Constants(
         tt, ss, ix, iy, iz, dxInv, dyInv, dzInv);
   }
}




/** Adapted from BlueBeats source code: FDLaplacian.h */
double
SaleheenDev::boundaryFDLaplacianSaleheen98SumPhi(const int index)
{
   double* phi = & (VmBlock_(index));
   const DiffusionCoefficients& diffConst = diffIntra_(index);
   
   double SumAphi =
        (diffConst.A1  * ( *(phi + xp1yz_)))
      + (diffConst.A2  * ( *(phi + xyp1z_)))
      + (diffConst.A3  * ( *(phi + xm1yz_)))
      + (diffConst.A4  * ( *(phi + xym1z_)))
      + (diffConst.A5  * ( *(phi + xp1yp1z_)))
      + (diffConst.A6  * ( *(phi + xm1yp1z_)))
      + (diffConst.A7  * ( *(phi + xm1ym1z_)))
      + (diffConst.A8  * ( *(phi + xp1ym1z_)))
      + (diffConst.A9  * ( *(phi + xyzp1_)))
      + (diffConst.A10 * ( *(phi + xyzm1_)))
      + (diffConst.A11 * ( *(phi + xyp1zp1_)))
      + (diffConst.A12 * ( *(phi + xyp1zm1_)))
      + (diffConst.A13 * ( *(phi + xym1zm1_)))
      + (diffConst.A14 * ( *(phi + xym1zp1_)))
      + (diffConst.A15 * ( *(phi + xp1yzp1_)))
      + (diffConst.A16 * ( *(phi + xm1yzp1_)))
      + (diffConst.A17 * ( *(phi + xm1yzm1_)))
      + (diffConst.A18 * ( *(phi + xp1yzm1_)));
  
   double result = SumAphi - (diffConst.sumA * (*phi));
   return result;
}

/** Make sure to use a one-sided difference approximation for the
 *  gradient of the conductivities for cells on the organ surface.
 *  (Conductivities are zero outside tissue cells and so cannot be used
 *  for a two-sided difference.)
 */
void
SaleheenDev::boundaryFDLaplacianSaleheen98Constants(
   const int*** tissue,
   const SymmetricTensor*** sigmaMatrix,
   const int& x, const int& y, const int& z,
   const double& dxInv, const double& dyInv, const double& dzInv)
{
   DiffusionCoefficients*** diffConst = diffIntra_.cArray();



   SymmetricTensor conductivityMatrix[3][3][3];
  
   for (int i=0; i<3; ++i)
      for (int j=0; j<3; ++j)
         for (int k=0; k<3; ++k)
         {
            conductivityMatrix[i][j][k].a11 = sigmaMatrix[x-1+i][y-1+j][z-1+k].a11;
            conductivityMatrix[i][j][k].a12 = sigmaMatrix[x-1+i][y-1+j][z-1+k].a12;
            conductivityMatrix[i][j][k].a13 = sigmaMatrix[x-1+i][y-1+j][z-1+k].a13;
            conductivityMatrix[i][j][k].a22 = sigmaMatrix[x-1+i][y-1+j][z-1+k].a22;
            conductivityMatrix[i][j][k].a23 = sigmaMatrix[x-1+i][y-1+j][z-1+k].a23;
            conductivityMatrix[i][j][k].a33 = sigmaMatrix[x-1+i][y-1+j][z-1+k].a33;
         }
                                                                                                                  
   const int i = 1;
   const int j = 1;
   const int k = 1;

   double dSigma11dx = dSigma(conductivityMatrix[i][j][k].a11, conductivityMatrix[i+1][j  ][k  ].a11, conductivityMatrix[i-1][j  ][k  ].a11, dxInv, tissue[x+i][y  ][z  ], tissue[x-1][y  ][z  ]);
   double dSigma12dy = dSigma(conductivityMatrix[i][j][k].a12, conductivityMatrix[i  ][j+1][k  ].a12, conductivityMatrix[i  ][j-1][k  ].a12, dyInv, tissue[x  ][y+1][z  ], tissue[x  ][y-1][z  ]);
   double dSigma13dz = dSigma(conductivityMatrix[i][j][k].a13, conductivityMatrix[i  ][j  ][k+1].a13, conductivityMatrix[i  ][j  ][k-1].a13, dzInv, tissue[x  ][y  ][z+1], tissue[x  ][y  ][z-1]);
   double dSigma12dx = dSigma(conductivityMatrix[i][j][k].a12, conductivityMatrix[i+1][j  ][k  ].a12, conductivityMatrix[i-1][j  ][k  ].a12, dxInv, tissue[x+1][y  ][z  ], tissue[x-1][y  ][z  ]);
   double dSigma22dy = dSigma(conductivityMatrix[i][j][k].a22, conductivityMatrix[i  ][j+1][k  ].a22, conductivityMatrix[i  ][j-1][k  ].a22, dyInv, tissue[x  ][y+1][z  ], tissue[x  ][y-1][z  ]);
   double dSigma23dz = dSigma(conductivityMatrix[i][j][k].a23, conductivityMatrix[i  ][j  ][k+1].a23, conductivityMatrix[i  ][j  ][k-1].a23, dzInv, tissue[x  ][y  ][z+1], tissue[x  ][y  ][z-1]);
   double dSigma13dx = dSigma(conductivityMatrix[i][j][k].a13, conductivityMatrix[i+1][j  ][k  ].a13, conductivityMatrix[i-1][j  ][k  ].a13, dxInv, tissue[x+1][y  ][z  ], tissue[x-1][y  ][z  ]);
   double dSigma23dy = dSigma(conductivityMatrix[i][j][k].a23, conductivityMatrix[i  ][j+1][k  ].a23, conductivityMatrix[i  ][j-1][k  ].a23, dyInv, tissue[x  ][y+1][z  ], tissue[x  ][y-1][z  ]);
   double dSigma33dz = dSigma(conductivityMatrix[i][j][k].a33, conductivityMatrix[i  ][j  ][k+1].a33, conductivityMatrix[i  ][j  ][k-1].a33, dzInv, tissue[x  ][y  ][z+1], tissue[x  ][y  ][z-1]);

   
   double sigmaX = 0.5 * dxInv * (dSigma11dx + dSigma12dy + dSigma13dz);
   double sigmaY = 0.5 * dyInv * (dSigma12dx + dSigma22dy + dSigma23dz);
   double sigmaZ = 0.5 * dzInv * (dSigma13dx + dSigma23dy + dSigma33dz);

   double dxdxInv = dxInv * dxInv;
   double dydyInv = dyInv * dyInv;
   double dzdzInv = dzInv * dzInv;
                      
   //ewd: for 1D or 2D case, check that a11,a22,a33 nonzero before divide:
   if (conductivityMatrix[i][j][k].a11 == 0.0)
   {
      diffConst[x][y][z].A1  = conductivityMatrix[i+1][j][k].a11 * dxdxInv;
      diffConst[x][y][z].A3  = conductivityMatrix[i-1][j][k].a11 * dxdxInv;
   }
   else
   {
      diffConst[x][y][z].A1  = conductivityMatrix[i+1][j][k].a11 * (dxdxInv - (sigmaX/conductivityMatrix[i][j][k].a11));
      diffConst[x][y][z].A3  = conductivityMatrix[i-1][j][k].a11 * (dxdxInv + (sigmaX/conductivityMatrix[i][j][k].a11));
   }
   if (conductivityMatrix[i][j][k].a22 == 0.0)
   {
      diffConst[x][y][z].A2  = conductivityMatrix[i][j+1][k].a22 * dydyInv;
      diffConst[x][y][z].A4  = conductivityMatrix[i][j-1][k].a22 * dydyInv;
   }
   else
   {
      diffConst[x][y][z].A2  = conductivityMatrix[i][j+1][k].a22 * (dydyInv - (sigmaY/conductivityMatrix[i][j][k].a22));
      diffConst[x][y][z].A4  = conductivityMatrix[i][j-1][k].a22 * (dydyInv + (sigmaY/conductivityMatrix[i][j][k].a22));
   }
   diffConst[x][y][z].A5  = diffConst[x][y][z].A6 = diffConst[x][y][z].A7 = diffConst[x][y][z].A8 = (0.5 * dxInv * dyInv);
  
   diffConst[x][y][z].A5  = conductivityMatrix[i+1][j+1][k].a12 * diffConst[x][y][z].A5;
   diffConst[x][y][z].A6  = -(conductivityMatrix[i-1][j+1][k].a12 * diffConst[x][y][z].A6);
   diffConst[x][y][z].A7  = conductivityMatrix[i-1][j-1][k].a12 * diffConst[x][y][z].A7;
   diffConst[x][y][z].A8  = -(conductivityMatrix[i+1][j-1][k].a12 * diffConst[x][y][z].A8);
   if (conductivityMatrix[i][j][k].a11 == 0.0)
   {
      diffConst[x][y][z].A9  = conductivityMatrix[i][j][k+1].a33 * dzdzInv;
      diffConst[x][y][z].A10 = conductivityMatrix[i][j][k-1].a33 * dzdzInv;
   }
   else
   {
      diffConst[x][y][z].A9  = conductivityMatrix[i][j][k+1].a33 * (dzdzInv - (sigmaZ/conductivityMatrix[i][j][k].a33));
      diffConst[x][y][z].A10 = conductivityMatrix[i][j][k-1].a33 * (dzdzInv + (sigmaZ/conductivityMatrix[i][j][k].a33));
   }
   diffConst[x][y][z].A11 = diffConst[x][y][z].A12 = diffConst[x][y][z].A13 = diffConst[x][y][z].A14 = (0.5 * dyInv * dzInv);
  
   diffConst[x][y][z].A11 = conductivityMatrix[i][j+1][k+1].a23 * diffConst[x][y][z].A11;
   diffConst[x][y][z].A12 = -(conductivityMatrix[i][j+1][k-1].a23 * diffConst[x][y][z].A12);
   diffConst[x][y][z].A13 = conductivityMatrix[i][j-1][k-1].a23 * diffConst[x][y][z].A13;
   diffConst[x][y][z].A14 = -(conductivityMatrix[i][j-1][k+1].a23 * diffConst[x][y][z].A14);

   diffConst[x][y][z].A15 = diffConst[x][y][z].A16 = diffConst[x][y][z].A17 = diffConst[x][y][z].A18 = (0.5 * dxInv * dzInv);

   diffConst[x][y][z].A15 = conductivityMatrix[i+1][j][k+1].a13 * diffConst[x][y][z].A15;
   diffConst[x][y][z].A16 = -(conductivityMatrix[i-1][j][k+1].a13 * diffConst[x][y][z].A16);
   diffConst[x][y][z].A17 = conductivityMatrix[i-1][j][k-1].a13 * diffConst[x][y][z].A17;
   diffConst[x][y][z].A18 = -(conductivityMatrix[i+1][j][k-1].a13 * diffConst[x][y][z].A18);

   diffConst[x][y][z].sumA = diffConst[x][y][z].A1 
      + diffConst[x][y][z].A2 
      + diffConst[x][y][z].A3 
      + diffConst[x][y][z].A4 
      + diffConst[x][y][z].A5 
      + diffConst[x][y][z].A6 
      + diffConst[x][y][z].A7 
      + diffConst[x][y][z].A8 
      + diffConst[x][y][z].A9 
      + diffConst[x][y][z].A10 
      + diffConst[x][y][z].A11 
      + diffConst[x][y][z].A12 
      + diffConst[x][y][z].A13 
      + diffConst[x][y][z].A14 
      + diffConst[x][y][z].A15 
      + diffConst[x][y][z].A16 
      + diffConst[x][y][z].A17 
      + diffConst[x][y][z].A18;
}


void SaleheenDev::updateVoltageBlock(
   const std::vector<double>& Vm)
{
   for (unsigned ii=0; ii<Vm.size(); ++ii)
   {
      int index = blockIndex_[ii];
      VmBlock_(index) = Vm[ii];
   }
}

