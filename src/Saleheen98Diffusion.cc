#include "Saleheen98Diffusion.hh"

#include "Anatomy.hh"
#include "SymmetricTensor.hh"
#include "DiffusionUtils.hh"
#include <cstdio> 
#include <iostream>

using namespace std;

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



Saleheen98PrecomputeDiffusion::Saleheen98PrecomputeDiffusion(
   const Saleheen98DiffusionParms& parms,
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
   
   diffIntra_.resize(nx, ny, nz);
   VmBlock_.resize(nx, ny, nz);

   buildTupleArray(anatomy);
   buildBlockIndex(anatomy);
   precomputeCoefficients(anatomy);
}

   
void Saleheen98PrecomputeDiffusion::calc(
   const vector<double>& Vm, vector<double>& dVm, double *recv_buf_, int nLocal)
{
   updateVoltageBlock(Vm, recv_buf_, nLocal);

   const int n = dVm.size();
#pragma omp parallel for
   for (int ii=0; ii<n; ++ii)
   {
      dVm[ii] = boundaryFDLaplacianSaleheen98SumPhi(localTuple_[ii]);
      dVm[ii] *= diffusionScale_;
   }
}

void Saleheen98PrecomputeDiffusion::calc_simd(
   const vector<double>& Vm, vector<double>& dVm, double *recv_buf_, int nLocal)
{
   updateVoltageBlock(Vm, recv_buf_, nLocal);
   uint32_t start = VmBlock_.tupleToIndex(1,1,1);
   uint32_t end = VmBlock_.tupleToIndex(localGrid_.nx()-1,localGrid_.ny(),localGrid_.nz());
   boundaryFDLaplacianSaleheen98SumPhi_All_simd(start,end-start,&(dVm[start]));
}

/** We're building the localTuple array only for local cells.  We can't
 * do stencil operations on remote particles so we shouldn't need
 * tuples.  We can use block indices instead.
 */
void Saleheen98PrecomputeDiffusion::buildTupleArray(const Anatomy& anatomy)
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

void Saleheen98PrecomputeDiffusion::buildBlockIndex(const Anatomy& anatomy)
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
Saleheen98PrecomputeDiffusion::precomputeCoefficients(const Anatomy& anatomy)
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

   //printAllConductivities(tissue, sigmaMintra); //ddt
   
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
  
   //printAllDiffusionWeights(tissue); //ddt
}




/** Adapted from BlueBeats source code: FDLaplacian.h */
double
Saleheen98PrecomputeDiffusion::boundaryFDLaplacianSaleheen98SumPhi(
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

// simdiazed version
// 'start' cannot be in the middle. 'start' must point to (n,0,0). 

void
Saleheen98PrecomputeDiffusion::boundaryFDLaplacianSaleheen98SumPhi_All_simd(const uint32_t start,const int32_t chunk_size, double* out)
{
  int ii;

  const unsigned Nx2 = localGrid_.nx();
  const unsigned Ny2 = localGrid_.ny();
  const unsigned Nz2 = localGrid_.nz();

  double* VmM = VmBlock_.cBlock();

  int xm1ym1z_ = ((-1) *Ny2 + (-1)) * Nz2;
  int xm1yz_ =   ((-1) *Ny2 + (0)) * Nz2;
  int xm1yp1z_ = ((-1) *Ny2 + (+1)) * Nz2;
  int xym1z_ =   ((0)  *Ny2 + (-1)) * Nz2;
  int xyp1z_ =   ((0)  *Ny2 + (+1)) * Nz2;
  int xp1ym1z_ = ((+1) *Ny2 + (-1)) * Nz2;
  int xp1yz_ =   ((+1) *Ny2 + (0)) * Nz2;
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
 
  simd_diff_ += start - 4;

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




/** Make sure to use a one-sided difference approximation for the
 *  gradient of the conductivities for cells on the organ surface.
 *  (Conductivities are zero outside tissue cells and so cannot be used
 *  for a two-sided difference.)
 */
void
Saleheen98PrecomputeDiffusion::boundaryFDLaplacianSaleheen98Constants(
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
   if (conductivityMatrix[i][j][k].a33 == 0.0)
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


void Saleheen98PrecomputeDiffusion::updateVoltageBlock(
   const vector<double>& Vm, double *recv_buf, int nLocal)
{
// hfwen
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

void Saleheen98PrecomputeDiffusion::printAllConductivities(
   const Array3d<int>& tissue, const Array3d<SymmetricTensor>& sigma)
{
   unsigned nx = localGrid_.nx();
   unsigned ny = localGrid_.ny();
   unsigned nz = localGrid_.nz();

   for (unsigned ix=0; ix<nx; ++ix)
      for (unsigned iy=0; iy<ny; ++iy)
         for (unsigned iz=0; iz<nz; ++iz)
         {
            Tuple globalTuple = localGrid_.globalTuple(Tuple(ix, iy, iz));
            printf("Conductivity: %5d %5d %5d %4d %18.12e %18.12e %18.12e %18.12e %18.12e %18.12e\n",
                   globalTuple.x(),
                   globalTuple.y(),
                   globalTuple.z(),
                   tissue(ix, iy, iz),
                   sigma(ix, iy, iz).a11,
                   sigma(ix, iy, iz).a22,
                   sigma(ix, iy, iz).a33,
                   sigma(ix, iy, iz).a12,
                   sigma(ix, iy, iz).a13,
                   sigma(ix, iy, iz).a23);
         }
}

void Saleheen98PrecomputeDiffusion::printAllDiffusionWeights(const Array3d<int>& tissue)
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
             diffIntra_(xx, yy, zz).A1,
             diffIntra_(xx, yy, zz).A2,
             diffIntra_(xx, yy, zz).A3,
             diffIntra_(xx, yy, zz).A4,
             diffIntra_(xx, yy, zz).A5,
             diffIntra_(xx, yy, zz).A6,
             diffIntra_(xx, yy, zz).A7,
             diffIntra_(xx, yy, zz).A8,
             diffIntra_(xx, yy, zz).A9,
             diffIntra_(xx, yy, zz).A10,
             diffIntra_(xx, yy, zz).A11,
             diffIntra_(xx, yy, zz).A12,
             diffIntra_(xx, yy, zz).A13,
             diffIntra_(xx, yy, zz).A14,
             diffIntra_(xx, yy, zz).A15,
             diffIntra_(xx, yy, zz).A16,
             diffIntra_(xx, yy, zz).A17,
             diffIntra_(xx, yy, zz).A18,
             diffIntra_(xx, yy, zz).sumA);
   }
}

void Saleheen98PrecomputeDiffusion::reorder_Coeff()
{
  uint32_t idx=0,ll,ii,jj,kk;
  const uint32_t Nx2 = localGrid_.nx();
  const uint32_t Ny2 = localGrid_.ny();
  const uint32_t Nz2 = localGrid_.nz();

  DiffusionCoefficients*** diffConst = diffIntra_.cArray();
  simd_diff_org_ = new double[Nx2*Ny2*Nz2+4]; //allocate but not aligend. 
  double* tgt = (double*)((uint64_t)simd_diff_org_ + 32 - ((uint64_t)simd_diff_org_)%32) ; //aligning
  simd_diff_ = tgt;

  idx--;
  for(ii=0;ii<Nx2;ii++)
  for(jj=0;jj<Ny2;jj++)
  for(kk=0;kk<Nz2;kk+=4)
  {
//    if(ii==1 && jj==1 && kk==0) printf("idx=%d\n",idx);
    idx++; for(ll=0;ll<4;ll++) tgt[idx*4+ll]=diffConst[ii][jj][(kk+ll-1+Nz2)%Nz2].A9;
    idx++; for(ll=0;ll<4;ll++) tgt[idx*4+ll]=diffConst[ii][jj][(kk+ll-1+Nz2)%Nz2].A14;
    idx++; for(ll=0;ll<4;ll++) tgt[idx*4+ll]=diffConst[ii][jj][(kk+ll-1+Nz2)%Nz2].A16;
    idx++; for(ll=0;ll<4;ll++) tgt[idx*4+ll]=diffConst[ii][jj][(kk+ll-1+Nz2)%Nz2].A11;
    idx++; for(ll=0;ll<4;ll++) tgt[idx*4+ll]=diffConst[ii][jj][(kk+ll-1+Nz2)%Nz2].A15;

    idx++; for(ll=0;ll<4;ll++) tgt[idx*4+ll]=diffConst[ii][jj][kk+ll].sumA * (-1);
    idx++; for(ll=0;ll<4;ll++) tgt[idx*4+ll]=diffConst[ii][jj][kk+ll].A8;
    idx++; for(ll=0;ll<4;ll++) tgt[idx*4+ll]=diffConst[ii][jj][kk+ll].A7;
    idx++; for(ll=0;ll<4;ll++) tgt[idx*4+ll]=diffConst[ii][jj][kk+ll].A6;
    idx++; for(ll=0;ll<4;ll++) tgt[idx*4+ll]=diffConst[ii][jj][kk+ll].A5;
    idx++; for(ll=0;ll<4;ll++) tgt[idx*4+ll]=diffConst[ii][jj][kk+ll].A4;
    idx++; for(ll=0;ll<4;ll++) tgt[idx*4+ll]=diffConst[ii][jj][kk+ll].A3;
    idx++; for(ll=0;ll<4;ll++) tgt[idx*4+ll]=diffConst[ii][jj][kk+ll].A2;
    idx++; for(ll=0;ll<4;ll++) tgt[idx*4+ll]=diffConst[ii][jj][kk+ll].A1;

    idx++; for(ll=0;ll<4;ll++) tgt[idx*4+ll]=diffConst[ii][jj][(kk+ll+1+Nz2)%Nz2].A10;
    idx++; for(ll=0;ll<4;ll++) tgt[idx*4+ll]=diffConst[ii][jj][(kk+ll+1+Nz2)%Nz2].A13;
    idx++; for(ll=0;ll<4;ll++) tgt[idx*4+ll]=diffConst[ii][jj][(kk+ll+1+Nz2)%Nz2].A17;
    idx++; for(ll=0;ll<4;ll++) tgt[idx*4+ll]=diffConst[ii][jj][(kk+ll+1+Nz2)%Nz2].A12;
    idx++; for(ll=0;ll<4;ll++) tgt[idx*4+ll]=diffConst[ii][jj][(kk+ll+1+Nz2)%Nz2].A18;
  }
}






