#include "Salheen98Diffusion.hh"

#include "Anatomy.hh"
#include "Conductivity.hh"
#include "conductivityFactory.hh"
#include <algorithm>
#include <iostream>

using namespace std;


// To Do
//
// 1.  What should be the default value of the cell type when the tissue
// block is created in precomputeCoefficients?


/** This implements (part of) the initialization conventions found at
 *  lines 1027 through 1076 of BlueBeats.cpp.  I.e., 9 everywhere,
 *  except 0 around the surface of the block.  (We'll overwrite this
 *  with actual cell types for the cells we have on this task later).
 *  This is a necessary initialization since the cells around the edge
 *  of the block may be completely synthetic.  For cells at the outer
 *  wall of the heart their may be no nbr cell defined in the anatomy,
 *  yet we need to have something there when we want to compute stencil
 *  operations. */
void initializeTissueBlock(Array3d<int>& tissue)
{
   for (unsigned ii=0; ii<tissue.size(); ++ii)
      tissue(ii) = 9;
   // cells on the outer face of the block are set to 0
   for (unsigned ii=0; ii<tissue.nx(); ++ii)
      for (unsigned jj=0; jj<tissue.ny(); ++jj)
	 for (unsigned kk=0; kk<tissue.nz(); ++kk)
	    if ( ii==0 || ii==tissue.nx()-1 ||
		 jj==0 || jj==tissue.nx()-1 ||
		 kk==0 || kk==tissue.nx()-1 )
	       tissue(ii, jj, kk) = 0;
}

/** We want to find the boundingBox such that any stencil point of any
 *  local atom is in the box.  It is not suffiient merely to iterate all
 *  of the local and remote atoms and find the maximum extent.  There
 *  may be local cells that are on the outer or inner walls of the
 *  heart.  Such cells will have no remote cells to satisfy their
 *  stencil.  Therefore, the safe bet is to iterate the local cells and
 *  add the stencil size in each direction.
 */
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


Salheen98PrecomputeDiffusion::Salheen98PrecomputeDiffusion(
   const Salheen98DiffusionParms& parms,
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
   conductivity_ = conductivityFactory(parms.conductivityName_);
   precomputeCoefficients(anatomy);
}

   
void Salheen98PrecomputeDiffusion::diffusion(
   const vector<double>& Vm, vector<double>& Istim)
{
   updateVoltageBlock(Vm);

   for (unsigned ii=0; ii<Istim.size(); ++ii)
   {
      Istim[ii] = boundaryFDLaplacianSaleheen98SumPhi(localTuple_[ii]);
      Istim[ii] *= diffusionScale_;
   }
}
/** We're building the localTuple array only for local cells.  We can't
 * do stencil operations on remote particles so we shouldn't need
 * tuples.  We can use block indices instead.
 */
void Salheen98PrecomputeDiffusion::buildTupleArray(const Anatomy& anatomy)
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

void Salheen98PrecomputeDiffusion::buildBlockIndex(const Anatomy& anatomy)
{
   blockIndex_.resize(anatomy.size());
   for (unsigned ii=0; ii<anatomy.size(); ++ii)
   {
      Tuple globalTuple = anatomy.globalTuple(ii);
      Tuple ll = localGrid_.localTuple(globalTuple);
      blockIndex_[ii] = VmBlock_.tupleToIndex(ll.x(), ll.y(), ll.z());
   }
}


/** BlueBeats precomputes the difusion coefficients as follows:
 *
 *  1.  The coefficients are stored in a 3D array of struct diffusion.
 *  the diffusion struct is defined in conductivity.h and contains 
 *  the 18 coefficients and a sum.
 * 
 *  2.  Memory for the 3D array is allocated by calling diffusion3D
 *  function with the size of the local brick as an argument on about
 *  line 753 of BlueBeats.cpp
 *
 *  BlueBeats actually uses 3 3D arrays:
 *  - tissue      (an array of ints containing the tissue type)
 *  - sigmaMintra (an array of struct sigmatensorMatrix)
 *  - diffIntra   (an array of struct diffusion)
 *
 *  3.  The conductivity matrix for each cell is calculated by calling
 *  CalcConductivityMatrixIBT.  This function is called once with theta
 *  and phi = 0 for every cell in the local block.  Later, as the
 *  orientation data is read from disk it is called again for every cell
 *  for which the angles have a defined value (i.e., not 255)
 *
 *  4.  The coefficients are computed by calling
 *  boundaryFDLaplacianSaleheen98Constants on about line 1989 of
 *  BlueBeats.cpp.  This calculation requires:
 *  - the tissue type matrix
 *  - the sigmatensorMatrix
 *  - reciprocal grid spacing
 *  Note that the range over which the calculation is performed needs a
 *  little bit of study.  It does not appear to be as simple as 0 to
 *  xMaxLocal (etc.)  This makes sense since you can only do diffusion
 *  for the local cells so there is no need to compute the coefficients
 *  for the edges of the local brick.  In fact, doing so will fail since
 *  the stencil can't be satisfied.
 
 */
void
Salheen98PrecomputeDiffusion::precomputeCoefficients(const Anatomy& anatomy)
{
   unsigned nx = localGrid_.nx();
   unsigned ny = localGrid_.ny();
   unsigned nz = localGrid_.nz();

   Array3d<SigmaTensorMatrix> sigmaMintra(nx, ny, nz);
   Array3d<int> tissue(nx, ny, nz);
   initializeTissueBlock(tissue);
   // What about default values for sigmaMintra and tissue?
   // Not all entries in the block are calculated.
   for (unsigned ii=0; ii<anatomy.size(); ++ii)
   {
      unsigned ib = blockIndex_[ii];
      int theta = anatomy.theta(ii);
      int phi = anatomy.phi(ii);
      conductivity_->compute(theta, phi, sigmaMintra(ib));
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
      const SigmaTensorMatrix*** ss =
	 (const SigmaTensorMatrix***) sigmaMintra.cArray();
      boundaryFDLaplacianSaleheen98Constants(
	 tt, ss, ix, iy, iz, dxInv, dyInv, dzInv);
   }
   
}




/** Adapted from BlueBeats source code: FDLaplacian.h */
double
Salheen98PrecomputeDiffusion::boundaryFDLaplacianSaleheen98SumPhi(
   const Tuple& tuple)
{
   DiffusionCoefficients*** diffConst = diffIntra_.cArray();
   double*** phi = VmBlock_.cArray();
   const unsigned& x = tuple.x();
   const unsigned& y = tuple.y();
   const unsigned& z = tuple.z();
   
   double SumAphi = (diffConst[x][y][z].A1  * (phi[x+1][y][z]))
                  + (diffConst[x][y][z].A2  * (phi[x][y+1][z]))
                  + (diffConst[x][y][z].A3  * (phi[x-1][y][z]))
                  + (diffConst[x][y][z].A4  * (phi[x][y-1][z]))
                  + (diffConst[x][y][z].A5  * (phi[x+1][y+1][z]))
                  + (diffConst[x][y][z].A6  * (phi[x-1][y+1][z]))
                  + (diffConst[x][y][z].A7  * (phi[x-1][y-1][z]))
                  + (diffConst[x][y][z].A8  * (phi[x+1][y-1][z]))
                  + (diffConst[x][y][z].A9  * (phi[x][y][z+1]))
                  + (diffConst[x][y][z].A10 * (phi[x][y][z-1]))
                  + (diffConst[x][y][z].A11 * (phi[x][y+1][z+1]))
                  + (diffConst[x][y][z].A12 * (phi[x][y+1][z-1]))
                  + (diffConst[x][y][z].A13 * (phi[x][y-1][z-1]))
                  + (diffConst[x][y][z].A14 * (phi[x][y-1][z+1]))
                  + (diffConst[x][y][z].A15 * (phi[x+1][y][z+1]))
                  + (diffConst[x][y][z].A16 * (phi[x-1][y][z+1]))
                  + (diffConst[x][y][z].A17 * (phi[x-1][y][z-1]))
                  + (diffConst[x][y][z].A18 * (phi[x+1][y][z-1]));
  
  double result = SumAphi - (diffConst[x][y][z].sumA * (phi[x][y][z]));
  return result;
}

/** Adapted from BlueBeats source code: FDLaplacian.h */
void
Salheen98PrecomputeDiffusion::boundaryFDLaplacianSaleheen98Constants(
   const int*** tissue,
   const SigmaTensorMatrix*** sigmaMatrix,
   const int& x, const int& y, const int& z,
   const double& dxInv, const double& dyInv, const double& dzInv)
{
   DiffusionCoefficients*** diffConst = diffIntra_.cArray();



   SigmaTensorMatrix conductivityMatrix[3][3][3];
  
   for (int i=0; i<3; ++i)
      for (int j=0; j<3; ++j)
	 for (int k=0; k<3; ++k)
	 {
	    conductivityMatrix[i][j][k].a11 = ((tissue[x-1+i][y-1+j][z-1+k]!=0)? (sigmaMatrix[x-1+i][y-1+j][z-1+k].a11):0.0);
	    conductivityMatrix[i][j][k].a12 = ((tissue[x-1+i][y-1+j][z-1+k]!=0)? (sigmaMatrix[x-1+i][y-1+j][z-1+k].a12):0.0);
	    conductivityMatrix[i][j][k].a13 = ((tissue[x-1+i][y-1+j][z-1+k]!=0)? (sigmaMatrix[x-1+i][y-1+j][z-1+k].a13):0.0);
	    conductivityMatrix[i][j][k].a22 = ((tissue[x-1+i][y-1+j][z-1+k]!=0)? (sigmaMatrix[x-1+i][y-1+j][z-1+k].a22):0.0);
	    conductivityMatrix[i][j][k].a23 = ((tissue[x-1+i][y-1+j][z-1+k]!=0)? (sigmaMatrix[x-1+i][y-1+j][z-1+k].a23):0.0);
	    conductivityMatrix[i][j][k].a33 = ((tissue[x-1+i][y-1+j][z-1+k]!=0)? (sigmaMatrix[x-1+i][y-1+j][z-1+k].a33):0.0);
	 }
                                                                                                                  
   double dxdxInv = dxInv * dxInv;
   double dydyInv = dyInv * dyInv;
   double dzdzInv = dzInv * dzInv;
                      
   const int i = 1;
   const int j = 1;
   const int k = 1;

   // original formulation - note that I replaced division by multiplication of reciprocal and division by 2 with multiplication by 0.5
   // I also used the central difference method for the differentiation. This results in a factor of 0.25 below
   // NOTE: sigmaX, sigmaY, sigmaZ could be precalculated if no deformation and no chance in conductivity
                                                                                                                    
   double sigmaX = 0.25 * dxInv * (((conductivityMatrix[i+1][j][k].a11 - conductivityMatrix[i-1][j][k].a11) * dxInv)
					+ ((conductivityMatrix[i][j+1][k].a12 - conductivityMatrix[i][j-1][k].a12) * dyInv)
					+ ((conductivityMatrix[i][j][k+1].a13 - conductivityMatrix[i][j][k-1].a13) * dzInv));
   double sigmaY = 0.25 * dyInv * (((conductivityMatrix[i+1][j][k].a12 - conductivityMatrix[i-1][j][k].a12) * dxInv)
					+ ((conductivityMatrix[i][j+1][k].a22 - conductivityMatrix[i][j-1][k].a12) * dyInv)
					+ ((conductivityMatrix[i][j][k+1].a23 - conductivityMatrix[i][j][k-1].a23) * dzInv));
   double sigmaZ = 0.25 * dzInv * (((conductivityMatrix[i+1][j][k].a13 - conductivityMatrix[i-1][j][k].a13) * dxInv)
					+ ((conductivityMatrix[i][j+1][k].a23 - conductivityMatrix[i][j-1][k].a23) * dyInv)
					+ ((conductivityMatrix[i][j][k+1].a33 - conductivityMatrix[i][j][k-1].a33) * dzInv));
  
   diffConst[x][y][z].A1  = conductivityMatrix[i+1][j][k].a11 * (dxdxInv - (sigmaX/conductivityMatrix[i][j][k].a11));
   diffConst[x][y][z].A3  = conductivityMatrix[i-1][j][k].a11 * (dxdxInv + (sigmaX/conductivityMatrix[i][j][k].a11));
   diffConst[x][y][z].A2  = conductivityMatrix[i][j+1][k].a22 * (dydyInv - (sigmaY/conductivityMatrix[i][j][k].a22));
   diffConst[x][y][z].A4  = conductivityMatrix[i][j-1][k].a22 * (dydyInv + (sigmaY/conductivityMatrix[i][j][k].a22));
   diffConst[x][y][z].A5  = diffConst[x][y][z].A6 = diffConst[x][y][z].A7 = diffConst[x][y][z].A8 = (0.5 * dxInv * dyInv);
  
   diffConst[x][y][z].A5  = conductivityMatrix[i+1][j+1][k].a12 * diffConst[x][y][z].A5;
   diffConst[x][y][z].A6  = -(conductivityMatrix[i-1][j+1][k].a12 * diffConst[x][y][z].A6);
   diffConst[x][y][z].A7  = conductivityMatrix[i-1][j-1][k].a12 * diffConst[x][y][z].A7;
   diffConst[x][y][z].A8  = -(conductivityMatrix[i+1][j-1][k].a12 * diffConst[x][y][z].A8);
   diffConst[x][y][z].A9  = conductivityMatrix[i][j][k+1].a33 * (dzdzInv - (sigmaZ/conductivityMatrix[i][j][k].a33));
   diffConst[x][y][z].A10 = conductivityMatrix[i][j][k-1].a33 * (dzdzInv + (sigmaZ/conductivityMatrix[i][j][k].a33));
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


void Salheen98PrecomputeDiffusion::updateVoltageBlock(
   const std::vector<double>& Vm)
{
   for (unsigned ii=0; ii<Vm.size(); ++ii)
   {
      int index = blockIndex_[ii];
      VmBlock_(index) = Vm[ii];
   }
}

